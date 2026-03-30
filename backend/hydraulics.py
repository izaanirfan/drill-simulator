import math
from trajectory import build_trajectory
from rheology import hb_fit
from temperature import temp_transient

def annular_area(Do_in, Di_in):
    return max((math.pi / 4) * (Do_in**2 - Di_in**2), 0.1)

def hydraulic_diameter(Do_in, Di_in):
    return max(Do_in - Di_in, 0.1)

def build_bha_profile(bha, depth):
    profile = []
    current_bottom = depth
    for comp in reversed(bha):
        length = float(getattr(comp, "length", 0))
        od = float(getattr(comp, "od", 5.0))
        top = max(current_bottom - length, 0)
        profile.append({"top": top, "bottom": current_bottom, "od": od})
        current_bottom = top
    return profile

def get_pipe_od(profile, md):
    for p in profile:
        if p["top"] <= md <= p["bottom"]:
            return p["od"]
    return 5.0  

def get_annulus_id(sections, md):
    for sec in sections:
        top_md = float(getattr(sec, "top_md", 0))
        end_md = float(getattr(sec, "end_md", 0))
        if top_md <= md <= end_md:
            t = str(getattr(sec, "type", "")).lower()
            return float(getattr(sec, "hole_d" if t == "open hole" else "casing_id", 8.5))
    return float(getattr(sections[-1], "hole_d", 8.5)) if sections else 8.5

def get_tvd(trajectory_profile, md):
    if not trajectory_profile: return md
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
        Q = float(getattr(data, "flowrate", 400))
        total_depth = max(float(getattr(data, "depth", 10000)), 100)
        
        fluid = getattr(data, "fluid", None)
        mw = float(getattr(fluid, "mw", 10)) if fluid else 10.0
        f600 = float(getattr(fluid, "fann_600", 60)) if fluid else 60.0
        f300 = float(getattr(fluid, "fann_300", 40)) if fluid else 40.0
        f200 = float(getattr(fluid, "fann_200", 30)) if fluid else 30.0
        f100 = float(getattr(fluid, "fann_100", 20)) if fluid else 20.0
        f6 = float(getattr(fluid, "fann_6", 10)) if fluid else 10.0
        f3 = float(getattr(fluid, "fann_3", 5)) if fluid else 5.0

        sbp = float(getattr(data, "sbp", 0.0)) 
        precision_const = 0.051948

        tau_y, K, n = hb_fit([f600, f300, f200, f100, f6, f3])
        n = max(min(n, 1.0), 0.1)
        K = max(K, 0.0001)
        
        trajectory_list = getattr(data, "trajectory", [])
        traj_profile = build_trajectory(trajectory_list)
        
        bha_list = getattr(data, "bha", [])
        bha_profile = build_bha_profile(bha_list, total_depth)
        
        well_sections = getattr(data, "well_sections", [])

        depths, ecd_profile, esd_profile, temp_profile = [], [], [], []
        cumulative_annular_loss = 0.0
        
        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))
        if md_points[-1] < total_depth: md_points.append(total_depth)

        temperature_obj = getattr(data, "temperature", None)
        Tsurf = float(getattr(temperature_obj, "surface_temp", 80)) if temperature_obj else 80.0
        Tbh = float(getattr(temperature_obj, "bhct", 150)) if temperature_obj else 150.0
        beta = float(getattr(temperature_obj, "beta", 0.0003)) if temperature_obj else 0.0003

        for md in md_points:
            Do = get_annulus_id(well_sections, md)
            Di = get_pipe_od(bha_profile, md)
            
            actual_step = md if len(depths) == 0 else md - depths[-1]
            dp = calculate_annular_pressure_loss(mw, Q, Do, Di, tau_y, K, n, actual_step)
            cumulative_annular_loss += dp

            tvd = get_tvd(traj_profile, md)
            if tvd <= 0.1: tvd = 0.1 

            ecd_base = mw + ((cumulative_annular_loss + sbp) / (precision_const * tvd))
            esd_base = mw + (sbp / (precision_const * tvd))

            T_current = temp_transient(md, Tsurf, Tbh, total_depth)
            temp_factor = 1 - beta * (T_current - Tsurf)

            depths.append(md)
            ecd_profile.append(round(ecd_base * temp_factor, 3))
            esd_profile.append(round(esd_base * temp_factor, 3))
            temp_profile.append(round(T_current, 1))

        # SPP Approximation
        p_bit = (mw * Q**2) / 10858 if Q > 0 else 0
        spp = p_bit + cumulative_annular_loss + (cumulative_annular_loss * 2.5)

        return {
            "summary": {
                "ecd_bottom": ecd_profile[-1] if ecd_profile else mw,
                "esd_bottom": esd_profile[-1] if esd_profile else mw,
                "total_annular_loss_psi": round(cumulative_annular_loss, 2),
                "spp": round(spp, 1)
            },
            "profile": {
                "depth": depths,
                "ecd": ecd_profile,
                "esd": esd_profile,
                "temp": temp_profile
            }
        }
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"error": str(e)}