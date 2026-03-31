import math

# Import your custom modules
from trajectory import build_trajectory
from rheology import hb_fit
from temperature import temp_transient

def get_val(obj, key, default=None):
    """Safely extracts data whether FastAPI passes a Pydantic object or a dictionary"""
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)

def build_bha_profile(bha_list, depth):
    profile = []
    current_bottom = depth
    for comp in reversed(bha_list):
        comp_len = float(get_val(comp, "length", 0))
        comp_od = float(get_val(comp, "od", 5.0))
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
        sec_top = float(get_val(sec, 'top_md', 0))
        sec_bot = float(get_val(sec, 'end_md', 0))
        sec_type = str(get_val(sec, 'type', '')).lower()
        if sec_top <= md <= sec_bot:
            if sec_type == "open hole":
                return float(get_val(sec, 'hole_d', 8.5))
            else:
                return float(get_val(sec, 'casing_id', 8.5))
    if sections:
        last_sec = sections[-1]
        return float(get_val(last_sec, 'hole_d', 8.5))
    return 8.5

def get_tvd(trajectory_profile, md):
    if not trajectory_profile:
        return md
    for i in range(len(trajectory_profile) - 1):
        p1 = trajectory_profile[i]
        p2 = trajectory_profile[i+1]
        if p1["md"] <= md <= p2["md"]:
            if p2["md"] == p1["md"]: return p1["tvd"]
            ratio = (md - p1["md"]) / (p2["md"] - p1["md"])
            return p1["tvd"] + ratio * (p2["tvd"] - p1["tvd"])
    return trajectory_profile[-1]["tvd"] if md >= trajectory_profile[-1]["md"] else md

def calculate_annular_pressure_loss(mw, flowrate, Do, Di, tau_y, K, n, length):
    if flowrate <= 0.1: return 0.0
    dh = max(Do - Di, 0.1)
    v_ft_min = (24.48 * flowrate) / (Do**2 - Di**2)
    v_ft_sec = v_ft_min / 60.0
    if v_ft_sec <= 0: return 0.0

    re_g = (89100 * mw * (v_ft_sec ** (2 - n))) / (K * ((144 / dh) ** n))
    re_g = max(re_g, 1.0) 
    re_critical = 3470 - (1370 * n)
    
    if re_g < re_critical:
        f = 24.0 / re_g 
    else:
        a = (math.log10(n) + 3.93) / 50.0
        b = (1.75 - math.log10(n)) / 7.0
        f = a / (re_g ** b)
    
    dp_psi = (f * length * mw * (v_ft_sec**2)) / (25.81 * dh)
    return max(dp_psi, 0.0)

def run_simulation(data):
    try:
        Q = float(get_val(data, 'flowrate', 400))
        total_depth = max(float(get_val(data, 'depth', 10000)), 100)
        sbp = float(get_val(data, 'sbp', 0.0))
        
        fluid = get_val(data, 'fluid', {})
        mw = float(get_val(fluid, 'mw', 10))
        
        fann_vals = [
            float(get_val(fluid, 'fann_600', 60)),
            float(get_val(fluid, 'fann_300', 40)),
            float(get_val(fluid, 'fann_200', 30)),
            float(get_val(fluid, 'fann_100', 20)),
            float(get_val(fluid, 'fann_6', 10)),
            float(get_val(fluid, 'fann_3', 5))
        ]
        tau_y, K, n = hb_fit(fann_vals)
        n = max(min(n, 1.0), 0.1)
        K = max(K, 0.0001)
        
        temp_data = get_val(data, 'temperature', {})
        t_surf = float(get_val(temp_data, 'surface_temp', 80))
        beta = float(get_val(temp_data, 'beta', 0.0003))
        geo_grad = float(get_val(temp_data, 'geothermal_grad', 1.5)) / 100.0
        
        t_rock_td = t_surf + (geo_grad * total_depth)
        cooling_factor = 0.85 if Q > 100 else 0.95
        t_bhct = t_surf + ((t_rock_td - t_surf) * cooling_factor)
        friction_heat_factor = 0.061 / max(mw, 1.0)
        
        traj_profile = build_trajectory(get_val(data, 'trajectory', []))
        bha_profile = build_bha_profile(get_val(data, 'bha', []), total_depth)
        well_sections = get_val(data, 'well_sections', [])

        depths, ecd_prof, esd_prof, temp_prof = [], [], [], []
        dp_segments = []
        
        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))
        if md_points[-1] < total_depth: md_points.append(total_depth)

        # PASS 1: Friction
        for md in md_points:
            Do = get_annulus_id(well_sections, md)
            Di = get_pipe_od(bha_profile, md)
            actual_step = md if len(depths) == 0 else md - depths[-1]
            dp = calculate_annular_pressure_loss(mw, Q, Do, Di, tau_y, K, n, actual_step)
            dp_segments.append(dp)
            depths.append(md)

        total_loss = sum(dp_segments)
        
        # PASS 2: Thermodynamics
        t_surf_return_base = t_surf + ((t_rock_td - t_surf) * 0.2)
        cumulative_dp = 0.0

        for i, md in enumerate(depths):
            cumulative_dp += dp_segments[i]
            
            norm_z = md / total_depth
            t_base_convective = t_surf_return_base + (t_bhct - t_surf_return_base) * (norm_z ** 0.6)
            friction_heat_added = (total_loss - cumulative_dp) * friction_heat_factor
            t_dyn_annulus = t_base_convective + friction_heat_added if Q > 0.1 else (t_surf + geo_grad * md)
            temp_prof.append(round(t_dyn_annulus, 1))

            temp_factor = 1 - beta * (t_dyn_annulus - t_surf)
            tvd = get_tvd(traj_profile, md)
            if tvd <= 0.1: tvd = 0.1

            ecd_base = mw + ((cumulative_dp + sbp) / (0.051948 * tvd))
            esd_base = mw + (sbp / (0.051948 * tvd))

            ecd_prof.append(round(ecd_base * temp_factor, 3))
            esd_prof.append(round(esd_base * temp_factor, 3))

        return {
            "summary": {
                "ecd_bottom": ecd_prof[-1] if ecd_prof else mw,
                "esd_bottom": esd_prof[-1] if esd_prof else mw,
                "total_annular_loss_psi": round(total_loss, 1)
            },
            "profile": {
                "depth": depths,
                "ecd": ecd_prof,
                "esd": esd_prof,
                "temp": temp_prof
            }
        }
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"error": str(e)}