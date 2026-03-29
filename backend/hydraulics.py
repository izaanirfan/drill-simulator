import math

# Import your custom modules
from trajectory import build_trajectory
from rheology import hb_fit
from temperature import temp_transient

# -----------------------------
# GEOMETRY
# -----------------------------
def annular_area(Do_in, Di_in):
    return max((math.pi / 4) * (Do_in**2 - Di_in**2), 0.1)

def hydraulic_diameter(Do_in, Di_in):
    return max(Do_in - Di_in, 0.1)

# -----------------------------
# BHA & WELL GEOMETRY
# -----------------------------
def build_bha_profile(bha_list, depth):
    profile = []
    current_bottom = depth
    for comp in reversed(bha_list):
        # Using dict access in case it's a standard parsed JSON dict
        comp_len = float(comp.get("length", 0)) if isinstance(comp, dict) else comp.length
        comp_od = float(comp.get("od", 5.0)) if isinstance(comp, dict) else comp.od
        
        top = max(current_bottom - comp_len, 0)
        profile.append({"top": top, "bottom": current_bottom, "od": comp_od})
        current_bottom = top
    return profile

def get_pipe_od(profile, md):
    for p in profile:
        if p["top"] <= md <= p["bottom"]:
            return p["od"]
    return 5.0  

def get_annulus_id(sections, md):
    for sec in sections:
        sec_top = float(sec.get("top_md", 0)) if isinstance(sec, dict) else sec.top_md
        sec_bot = float(sec.get("end_md", 0)) if isinstance(sec, dict) else sec.end_md
        sec_type = sec.get("type", "").lower() if isinstance(sec, dict) else sec.type.lower()
        
        if sec_top <= md <= sec_bot:
            if sec_type == "open hole":
                return float(sec.get("hole_d", 8.5) if isinstance(sec, dict) else sec.hole_d)
            else:
                return float(sec.get("casing_id", 8.5) if isinstance(sec, dict) else sec.casing_id)
                
    # Fallback to last section
    last_sec = sections[-1]
    return float(last_sec.get("hole_d", 8.5) if isinstance(last_sec, dict) else last_sec.hole_d)

# -----------------------------
# TVD INTERPOLATION
# -----------------------------
def get_tvd(trajectory_profile, md):
    if not trajectory_profile:
        return md
    
    # Interpolate TVD between survey points
    for i in range(len(trajectory_profile) - 1):
        p1 = trajectory_profile[i]
        p2 = trajectory_profile[i+1]
        if p1["md"] <= md <= p2["md"]:
            if p2["md"] == p1["md"]:
                return p1["tvd"]
            ratio = (md - p1["md"]) / (p2["md"] - p1["md"])
            return p1["tvd"] + ratio * (p2["tvd"] - p1["tvd"])
            
    # Fallback if MD is past the last survey point
    return trajectory_profile[-1]["tvd"] if md >= trajectory_profile[-1]["md"] else md

# -----------------------------
# TRUE HERSCHEL-BULKLEY PRESSURE LOSS MODEL
# -----------------------------
def calculate_annular_pressure_loss(mw, flowrate, Do, Di, tau_y, K, n, length):
    """
    Calculates frictional pressure drop using the Highly Accurate 
    Generalized Reynolds Number (Re_G) for Yield-Power Law fluids.
    """
    # ZERO FLOW BYPASS: Prevents mathematical crashes for static EMW calculations
    if flowrate <= 0.1:
        return 0.0
        
    dh = max(Do - Di, 0.1)
    
    # Annular velocity in ft/sec
    v_ft_min = (24.48 * flowrate) / (Do**2 - Di**2)
    v_ft_sec = v_ft_min / 60.0
    
    if v_ft_sec <= 0:
        return 0.0

    # Generalized Reynolds Number (Re_G) for Non-Newtonian Fluids
    # 89100 is the standard API constant for oilfield units
    re_g = (89100 * mw * (v_ft_sec ** (2 - n))) / (K * ((144 / dh) ** n))
    re_g = max(re_g, 1.0) # Prevent zero division
    
    # Calculate Critical Reynolds Number (Point of turbulent transition)
    re_critical = 3470 - (1370 * n)
    
    # Friction Factor Selection
    if re_g < re_critical:
        # LAMINAR FLOW REGIME
        f = 24.0 / re_g 
    else:
        # TURBULENT FLOW REGIME (Dodge-Metzner empirical correlation)
        a = (math.log10(n) + 3.93) / 50.0
        b = (1.75 - math.log10(n)) / 7.0
        f = a / (re_g ** b)
    
    # Final Pressure Drop for this segment
    dp_psi = (f * length * mw * (v_ft_sec**2)) / (25.81 * dh)
    
    return max(dp_psi, 0.0)

# -----------------------------
# MAIN SIMULATION
# -----------------------------
def run_simulation(data):
    try:
        # Standardized dictionary access to prevent 500 errors from missing JSON fields
        Q = float(data.get('flowrate', 400))
        total_depth = max(float(data.get('depth', 10000)), 100)
        sbp = float(data.get('sbp', 0.0))
        
        fluid = data.get('fluid', {})
        mw = float(fluid.get('mw', 10))
        f600 = float(fluid.get('fann_600', 60))
        f300 = float(fluid.get('fann_300', 40))
        f200 = float(fluid.get('fann_200', 30))
        f100 = float(fluid.get('fann_100', 20))
        f6 = float(fluid.get('fann_6', 10))
        f3 = float(fluid.get('fann_3', 5))
        
        precision_const = 0.051948

        # 1. Use external rheology.py for HB Fit
        fann_readings = [f600, f300, f200, f100, f6, f3]
        tau_y, K, n = hb_fit(fann_readings)
        
        # Protect fluid parameters from extreme math bounds
        n = max(min(n, 1.0), 0.1)
        K = max(K, 0.0001)
        
        # 2. Build trajectory using trajectory.py
        traj_profile = build_trajectory(data.get('trajectory', []))
        bha_profile = build_bha_profile(data.get('bha', []), total_depth)

        depths, ecd_profile, esd_profile = [], [], []
        cumulative_annular_loss = 0.0
        
        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))

        for md in md_points:
            if md > total_depth: md = total_depth

            Do = get_annulus_id(data.get('well_sections', []), md)
            Di = get_pipe_od(bha_profile, md)
            
            # Calculate friction using the UPGRADED Herschel-Bulkley model
            dp = calculate_annular_pressure_loss(mw, Q, Do, Di, tau_y, K, n, step)
            cumulative_annular_loss += dp

            # 3. Interpolate EXACT TVD for this step
            tvd = get_tvd(traj_profile, md)
            if tvd <= 0.1: 
                tvd = 0.1 # Prevent Division by Zero at surface

            # 4. ECD/ESD = MW + ((Sum of dP + SBP) / (Constant * TVD))
            ecd_base = mw + ((cumulative_annular_loss + sbp) / (precision_const * tvd))
            esd_base = mw + (sbp / (precision_const * tvd))

            # 5. Use temperature.py for transient temp correction
            temp_data = data.get('temperature', {})
            Tsurf = float(temp_data.get('surface_temp', 80))
            Tbh = float(temp_data.get('bhct', 150))
            beta = float(temp_data.get('beta', 0.0003))
            
            T_current = temp_transient(md, Tsurf, Tbh, total_depth)
            temp_factor = 1 - beta * (T_current - Tsurf)

            depths.append(md)
            ecd_profile.append(round(ecd_base * temp_factor, 3))
            esd_profile.append(round(esd_base * temp_factor, 3))

        return {
            "summary": {
                "ecd_bottom": ecd_profile[-1] if ecd_profile else mw,
                "esd_bottom": esd_profile[-1] if esd_profile else mw,
                "total_annular_loss_psi": round(cumulative_annular_loss, 2)
            },
            "profile": {
                "depth": depths,
                "ecd": ecd_profile,
                "esd": esd_profile
            }
        }
    except Exception as e:
        return {"error": str(e)}