import math

# Import your custom modules
from trajectory import build_trajectory
from rheology import hb_fit

# -----------------------------
# HELPER: Data Extractor
# -----------------------------
def get_val(obj, key, default=None):
    """Safely extracts data whether FastAPI passes a Pydantic object or a dictionary"""
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)

# -----------------------------
# GEOMETRY
# -----------------------------
def annular_area(Do_in, Di_in):
    return max((math.pi / 4) * (Do_in**2 - Di_in**2), 0.1)

def hydraulic_diameter(Do_in, Di_in):
    return max(Do_in - Di_in, 0.1)

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

# -----------------------------
# HYDRAULICS ENGINE
# -----------------------------
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

# SINGLE PASS SOLVER
def run_pass(Q, sbp, total_depth, md_points, data, mw, tau_y, K, n, traj_profile, bha_profile, t_surf, t_bhct, t_rock_td, friction_heat_factor, beta, geo_grad):
    depths, ecd_prof, esd_prof, temp_prof = [], [], [], []
    cumulative_dp = 0.0
    
    # 1. Friction Pass
    dp_segments = []
    well_sections = get_val(data, 'well_sections', [])
    for md in md_points:
        Do = get_annulus_id(well_sections, md)
        Di = get_pipe_od(bha_profile, md)
        actual_step = md if len(depths) == 0 else md - depths[-1]
        dp = calculate_annular_pressure_loss(mw, Q, Do, Di, tau_y, K, n, actual_step)
        dp_segments.append(dp)
        depths.append(md)

    total_ann_loss = sum(dp_segments)
    
    # 2. Thermo Pass
    t_surf_return_base = t_surf + ((t_rock_td - t_surf) * 0.2)
    
    for i, md in enumerate(depths):
        cumulative_dp += dp_segments[i]
        
        norm_z = md / total_depth
        t_base_convective = t_surf_return_base + (t_bhct - t_surf_return_base) * (norm_z ** 0.6)
        friction_heat_added = (total_ann_loss - cumulative_dp) * friction_heat_factor
        t_dyn_annulus = t_base_convective + friction_heat_added if Q > 0.1 else (t_surf + geo_grad * md)
        temp_prof.append(round(t_dyn_annulus, 1))

        temp_factor = 1 - beta * (t_dyn_annulus - t_surf)
        tvd = get_tvd(traj_profile, md)
        if tvd <= 0.1: tvd = 0.1

        ecd_base = mw + ((cumulative_dp + sbp) / (0.051948 * tvd))
        esd_base = mw + (sbp / (0.051948 * tvd))

        ecd_prof.append(round(ecd_base * temp_factor, 3))
        esd_prof.append(round(esd_base * temp_factor, 3))
        
    return depths, ecd_prof, esd_prof, temp_prof, total_ann_loss

# -----------------------------
# MAIN SIMULATION EXPORT
# -----------------------------
def run_simulation(data):
    try:
        Q = float(get_val(data, 'flowrate', 400))
        total_depth = max(float(get_val(data, 'depth', 10000)), 100)
        
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

        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))
        if md_points[-1] < total_depth: md_points.append(total_depth)

        # -----------------------------
        # PROPRIETARY TWO-PASS SOLVER
        # -----------------------------
        main_mode = get_val(data, 'main_mode', 'conventional')
        mpd_mode = get_val(data, 'mpd_mode', 'sbp')
        applied_sbp = 0.0
        target_warning = None

        if main_mode == 'mpd' and mpd_mode == 'downhole':
            # PASS 1: Simulate at 0 SBP to find friction
            depths, ecd0, esd0, tmp0, loss0 = run_pass(Q, 0.0, total_depth, md_points, data, mw, tau_y, K, n, traj_profile, bha_profile, t_surf, t_bhct, t_rock_td, friction_heat_factor, beta, geo_grad)
            
            anchor_type = get_val(data, 'anchor_type', 'bh')
            if anchor_type == 'shoe':
                secs = get_val(data, 'well_sections', [])
                def get_type(s): return get_val(s, 'type', '').lower()
                def get_top(s): return float(get_val(s, 'top_md', total_depth))
                anchor_depth = next((get_top(s) for s in secs if get_type(s) == 'open hole'), total_depth)
            elif anchor_type == 'fixed':
                anchor_depth = float(get_val(data, 'anchor_fixed_val', 5000))
            else:
                anchor_depth = total_depth
                
            idx = next((i for i, d in enumerate(depths) if d >= anchor_depth), -1)
            base_ecd = ecd0[idx] if idx != -1 else ecd0[-1]
            
            target_type = get_val(data, 'target_type', 'density')
            target_val = float(get_val(data, 'target_val', 12.0))
            targ_press = (target_val * anchor_depth * 0.051948) if target_type == 'density' else target_val
            base_press = base_ecd * anchor_depth * 0.051948
            req_sbp = targ_press - base_press
            
            min_sbp = float(get_val(data, 'min_sbp', 0))
            max_sbp = float(get_val(data, 'max_sbp', 2000))
            
            if req_sbp < min_sbp:
                applied_sbp = min_sbp
                target_warning = f"Target EMW cannot be achieved (Requires SBP < Min SBP). Simulating using Min SBP ({min_sbp})."
            elif req_sbp > max_sbp:
                applied_sbp = max_sbp
                target_warning = f"Target EMW cannot be achieved (Requires SBP > Max SBP). Simulating using Max SBP ({max_sbp})."
            else:
                applied_sbp = req_sbp
        elif main_mode == 'mpd' and mpd_mode == 'sbp':
            applied_sbp = float(get_val(data, 'sbp', 0.0))

        # PASS 2: Final Run
        depths, ecd_prof, esd_prof, temp_prof, total_loss = run_pass(Q, applied_sbp, total_depth, md_points, data, mw, tau_y, K, n, traj_profile, bha_profile, t_surf, t_bhct, t_rock_td, friction_heat_factor, beta, geo_grad)

        # -----------------------------
        # STANDPIPE PRESSURE CALC
        # -----------------------------
        bit_tfa = float(get_val(data, 'bit_tfa', 1.0))
        p_bit = (mw * Q**2) / (10858 * bit_tfa**2) if Q > 0 else 0
        p_pipe = total_loss * 2.5 # Internal friction est
        spp = p_bit + total_loss + p_pipe

        # -----------------------------
        # ZERO FLOW THERMAL CALC
        # -----------------------------
        zf_mode = get_val(data, 'zero_flow_mode', 'immediate')
        zf_time = float(get_val(data, 'zero_flow_time', 12)) if zf_mode == 'insitu' else 0
        thermal_factor = 1 - math.exp(-zf_time / 10.0)
        
        cbhp_dyn = ecd_prof[-1] * 0.051948 * total_depth
        t_avg_geo_td = (t_surf + t_rock_td) / 2
        esd_geo_td = mw * (1 - beta * (t_avg_geo_td - t_surf))
        rho_t_td = esd_prof[-1] + (esd_geo_td - esd_prof[-1]) * thermal_factor
        
        zf_req_sbp = max(cbhp_dyn - (rho_t_td * 0.051948 * total_depth), 0)
        
        zf_emw_prof, t_geo_prof = [], []
        for i, md in enumerate(depths):
            t_geo = t_surf + geo_grad * md
            t_geo_prof.append(round(t_geo, 1))
            if md < 0.1:
                zf_emw_prof.append(mw)
            else:
                t_avg = (t_surf + t_geo) / 2
                esd_geo = mw * (1 - beta * (t_avg - t_surf))
                rho_t = esd_prof[i] + (esd_geo - esd_prof[i]) * thermal_factor
                hydro = rho_t * 0.051948 * md
                zf_emw_prof.append(round((zf_req_sbp + hydro) / (0.051948 * md), 3))

        return {
            "summary": {
                "ecd_bottom": ecd_prof[-1] if ecd_prof else mw,
                "esd_bottom": esd_prof[-1] if esd_prof else mw,
                "total_annular_loss_psi": round(total_loss, 1),
                "applied_sbp": round(applied_sbp, 1),
                "spp": round(spp, 1),
                "zf_req_sbp": round(zf_req_sbp, 1),
                "target_warning": target_warning
            },
            "profile": {
                "depth": depths,
                "ecd": ecd_prof,
                "esd": esd_prof,
                "temp": temp_prof,
                "zf_emw": zf_emw_prof,
                "t_geo": t_geo_prof
            }
        }
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"error": str(e)}