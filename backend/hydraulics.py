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
def build_bha_profile(bha, depth):
    profile = []
    current_bottom = depth
    for comp in reversed(bha):
        top = max(current_bottom - comp.length, 0)
        profile.append({"top": top, "bottom": current_bottom, "od": comp.od})
        current_bottom = top
    return profile

def get_pipe_od(profile, md):
    for p in profile:
        if p["top"] <= md <= p["bottom"]:
            return p["od"]
    return 5.0  

def get_annulus_id(sections, md):
    for sec in sections:
        if sec.top_md <= md <= sec.end_md:
            t = sec.type.lower()
            return sec.hole_d if t == "open hole" else sec.casing_id
    return sections[-1].hole_d

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
# REFINED PRESSURE LOSS MODEL
# -----------------------------
def calculate_annular_pressure_loss(mw, flowrate, Do, Di, mu_cp, length):
    dh = Do - Di
    v = (24.48 * flowrate) / (Do**2 - Di**2)
    reynolds = (15.47 * mw * v * dh) / mu_cp
    
    if reynolds < 2100:
        f = 64 / reynolds 
    else:
        f = 0.3164 / (reynolds**0.25) 
    
    dp_field = (f * length * mw * v**2) / (92000 * dh)
    return max(dp_field, 0)

# -----------------------------
# MAIN SIMULATION
# -----------------------------
def run_simulation(data):
    try:
        Q = float(data.flowrate or 400)
        total_depth = max(float(data.depth or 10000), 100)
        mw = float(data.fluid.mw or 10)
        sbp = float(data.sbp or 0.0) 
        precision_const = 0.051948

        # 1. Use external rheology.py for HB Fit
        fann_readings = [
            data.fluid.fann_600, data.fluid.fann_300, data.fluid.fann_200, 
            data.fluid.fann_100, data.fluid.fann_6, data.fluid.fann_3
        ]
        tau_y, K, n = hb_fit(fann_readings)
        
        # 2. Build trajectory using trajectory.py
        traj_profile = build_trajectory(data.trajectory)
        bha_profile = build_bha_profile(data.bha, total_depth)

        depths, ecd_profile, esd_profile = [], [], []
        cumulative_annular_loss = 0
        
        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))

        for md in md_points:
            if md > total_depth: md = total_depth

            Do = get_annulus_id(data.well_sections, md)
            Di = get_pipe_od(bha_profile, md)
            dh = hydraulic_diameter(Do, Di)
            
            v_ftmin = (24.48 * Q) / (Do**2 - Di**2)
            gamma = (1.6 * v_ftmin) / dh 
            mu_cp = ((tau_y + K * (gamma**n)) / gamma) * 478.8 if gamma > 0.1 else 100
            
            dp = calculate_annular_pressure_loss(mw, Q, Do, Di, mu_cp, step)
            cumulative_annular_loss += dp

            # 3. Interpolate EXACT TVD for this step
            tvd = get_tvd(traj_profile, md)
            if tvd <= 0.1: 
                tvd = 0.1 # Prevent Division by Zero at surface

            # 4. ECD/ESD = MW + ((Sum of dP + SBP) / (Constant * TVD))
            ecd_base = mw + ((cumulative_annular_loss + sbp) / (precision_const * tvd))
            esd_base = mw + (sbp / (precision_const * tvd))

            # 5. Use temperature.py for transient temp correction
            Tsurf = data.temperature.surface_temp
            Tbh = data.temperature.bhct
            T_current = temp_transient(md, Tsurf, Tbh, total_depth)
            temp_factor = 1 - data.temperature.beta * (T_current - Tsurf)

            depths.append(md)
            ecd_profile.append(round(ecd_base * temp_factor, 3))
            esd_profile.append(round(esd_base * temp_factor, 3))

        return {
            "summary": {
                "ecd_bottom": ecd_profile[-1],
                "esd_bottom": esd_profile[-1],
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