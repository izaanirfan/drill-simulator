import math
from trajectory import build_trajectory
from rheology import hb_fit

def annular_area(Do_in, Di_in):
    return max((math.pi / 4) * (Do_in**2 - Di_in**2), 0.1)

def hydraulic_diameter(Do_in, Di_in):
    return max(Do_in - Di_in, 0.1)

def build_bha_profile(bha, depth):
    profile = []
    current_bottom = depth
    for comp in reversed(bha):
        top = max(current_bottom - float(comp.length), 0)
        profile.append({"top": top, "bottom": current_bottom, "od": float(comp.od)})
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
            return float(sec.hole_d if t == "open hole" else sec.casing_id)
    return float(sections[-1].hole_d) if sections else 8.5

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
        Q = float(data.flowrate or 400)
        total_depth = max(float(data.depth or 10000), 100)
        mw = float(data.fluid.mw or 10)
        sbp = float(data.sbp or 0.0) 
        precision_const = 0.051948

        fann_readings = [
            data.fluid.fann_600, data.fluid.fann_300, data.fluid.fann_200, 
            data.fluid.fann_100, data.fluid.fann_6, data.fluid.fann_3
        ]
        tau_y, K, n = hb_fit(fann_readings)
        n = max(min(n, 1.0), 0.1)
        K = max(K, 0.0001)
        
        traj_profile = build_trajectory(data.trajectory)
        bha_profile = build_bha_profile(data.bha, total_depth)

        depths, ecd_profile, esd_profile, temp_profile = [], [], [], []
        cumulative_annular_loss = 0.0
        
        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))
        if md_points[-1] < total_depth:
            md_points.append(total_depth)

        # Thermal Defaults
        Tsurf = float(data.temperature.surface_temp or 80)
        Tbh = float(data.temperature.bhct or 150)
        beta = getattr(data.temperature, 'beta', 0.0003)
        geo_grad = (Tbh - Tsurf) / max(total_depth, 1)

        t_rock_td = Tsurf + (geo_grad * total_depth)
        cooling_factor = 0.85 if Q > 100 else 0.95
        t_bhct_dyn = Tsurf + ((t_rock_td - Tsurf) * cooling_factor)
        
        friction_heat_factor = 0.061 / max(mw, 1.0)
        t_surf_return_base = Tsurf + ((t_rock_td - Tsurf) * 0.2)

        for md in md_points:
            Do = get_annulus_id(data.well_sections, md)
            Di = get_pipe_od(bha_profile, md)
            
            actual_step = md if len(depths) == 0 else md - depths[-1]
            dp = calculate_annular_pressure_loss(mw, Q, Do, Di, tau_y, K, n, actual_step)
            cumulative_annular_loss += dp

            tvd = get_tvd(traj_profile, md)
            if tvd <= 0.1: tvd = 0.1 

            ecd_base = mw + ((cumulative_annular_loss + sbp) / (precision_const * tvd))
            esd_base = mw + (sbp / (precision_const * tvd))

            norm_z = md / total_depth
            t_base_convective = t_surf_return_base + (t_bhct_dyn - t_surf_return_base) * (norm_z ** 0.6)
            
            t_dyn_annulus = t_base_convective + (cumulative_annular_loss * friction_heat_factor) if Q > 0.1 else (Tsurf + geo_grad * md)

            temp_factor = 1 - beta * (t_dyn_annulus - Tsurf)

            depths.append(md)
            ecd_profile.append(round(ecd_base * temp_factor, 3))
            esd_profile.append(round(esd_base * temp_factor, 3))
            temp_profile.append(round(t_dyn_annulus, 1))

        # SPP Calc (Approximated if not passed)
        p_bit = (mw * Q**2) / 10858 if Q > 0 else 0
        p_pipe = cumulative_annular_loss * 2.5 
        spp = p_bit + cumulative_annular_loss + p_pipe

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