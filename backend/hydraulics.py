from flask import Flask, request, jsonify, send_from_directory
import math
import os

app = Flask(__name__, static_folder='static')

@app.route('/')
def index():
    return send_from_directory(app.static_folder, 'index.html')

# Import your custom modules
from trajectory import build_trajectory
from rheology import hb_fit

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
        comp_len = float(comp.get("length", 0)) if isinstance(comp, dict) else float(comp.length)
        comp_od = float(comp.get("od", 5.0)) if isinstance(comp, dict) else float(comp.od)
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
        sec_top = float(sec.get("top_md", 0)) if isinstance(sec, dict) else float(sec.top_md)
        sec_bot = float(sec.get("end_md", 0)) if isinstance(sec, dict) else float(sec.end_md)
        sec_type = sec.get("type", "").lower() if isinstance(sec, dict) else sec.type.lower()
        if sec_top <= md <= sec_bot:
            if sec_type == "open hole":
                return float(sec.get("hole_d", 8.5) if isinstance(sec, dict) else sec.hole_d)
            else:
                return float(sec.get("casing_id", 8.5) if isinstance(sec, dict) else sec.casing_id)
    last_sec = sections[-1] if sections else {}
    return float(last_sec.get("hole_d", 8.5) if isinstance(last_sec, dict) else 8.5)

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
    for md in md_points:
        Do = get_annulus_id(data.get('well_sections', []), md)
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
# MAIN SIMULATION API
# -----------------------------
@app.route('/simulate', methods=['POST'])
def run_simulation():
    try:
        data = request.json
        Q = float(data.get('flowrate', 400))
        total_depth = max(float(data.get('depth', 10000)), 100)
        
        fluid = data.get('fluid', {})
        mw = float(fluid.get('mw', 10))
        tau_y, K, n = hb_fit([float(fluid.get(f'fann_{x}', 0)) for x in [600, 300, 200, 100, 6, 3]])
        n = max(min(n, 1.0), 0.1)
        K = max(K, 0.0001)
        
        temp_data = data.get('temperature', {})
        t_surf = float(temp_data.get('surface_temp', 80))
        beta = float(temp_data.get('beta', 0.0003))
        geo_grad = float(temp_data.get('geothermal_grad', 1.5)) / 100.0
        
        t_rock_td = t_surf + (geo_grad * total_depth)
        cooling_factor = 0.85 if Q > 100 else 0.95
        t_bhct = t_surf + ((t_rock_td - t_surf) * cooling_factor)
        friction_heat_factor = 0.061 / max(mw, 1.0)
        
        traj_profile = build_trajectory(data.get('trajectory', []))
        bha_profile = build_bha_profile(data.get('bha', []), total_depth)

        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))
        if md_points[-1] < total_depth: md_points.append(total_depth)

        # -----------------------------
        # PROPRIETARY TWO-PASS SOLVER
        # -----------------------------
        main_mode = data.get('main_mode', 'conventional')
        mpd_mode = data.get('mpd_mode', 'sbp')
        applied_sbp = 0.0
        target_warning = None

        if main_mode == 'mpd' and mpd_mode == 'downhole':
            # PASS 1: Simulate at 0 SBP to find friction
            depths, ecd0, esd0, tmp0, loss0 = run_pass(Q, 0.0, total_depth, md_points, data, mw, tau_y, K, n, traj_profile, bha_profile, t_surf, t_bhct, t_rock_td, friction_heat_factor, beta, geo_grad)
            
            # Find Anchor Depth
            anchor_type = data.get('anchor_type', 'bh')
            if anchor_type == 'shoe':
                secs = data.get('well_sections', [])
                anchor_depth = next((float(s['top_md']) for s in secs if s['type'].lower() == 'open hole'), total_depth)
            elif anchor_type == 'fixed':
                anchor_depth = float(data.get('anchor_fixed_val', 5000))
            else:
                anchor_depth = total_depth
                
            # Interpolate ECD at Anchor
            idx = next((i for i, d in enumerate(depths) if d >= anchor_depth), -1)
            base_ecd = ecd0[idx] if idx != -1 else ecd0[-1]
            
            # Calculate required SBP
            target_type = data.get('target_type', 'density')
            target_val = float(data.get('target_val', 12.0))
            targ_press = (target_val * anchor_depth * 0.051948) if target_type == 'density' else target_val
            base_press = base_ecd * anchor_depth * 0.051948
            req_sbp = targ_press - base_press
            
            # Apply Limits
            min_sbp = float(data.get('min_sbp', 0))
            max_sbp = float(data.get('max_sbp', 2000))
            
            if req_sbp < min_sbp:
                applied_sbp = min_sbp
                target_warning = f"Target EMW cannot be achieved (Requires SBP < Min SBP). Simulating using Min SBP ({min_sbp})."
            elif req_sbp > max_sbp:
                applied_sbp = max_sbp
                target_warning = f"Target EMW cannot be achieved (Requires SBP > Max SBP). Simulating using Max SBP ({max_sbp})."
            else:
                applied_sbp = req_sbp
        elif main_mode == 'mpd' and mpd_mode == 'sbp':
            applied_sbp = float(data.get('sbp', 0.0))

        # PASS 2: Final Run
        depths, ecd_prof, esd_prof, temp_prof, total_loss = run_pass(Q, applied_sbp, total_depth, md_points, data, mw, tau_y, K, n, traj_profile, bha_profile, t_surf, t_bhct, t_rock_td, friction_heat_factor, beta, geo_grad)

        # -----------------------------
        # STANDPIPE PRESSURE CALC
        # -----------------------------
        bit_tfa = float(data.get('bit_tfa', 1.0))
        p_bit = (mw * Q**2) / (10858 * bit_tfa**2) if Q > 0 else 0
        p_pipe = total_loss * 2.5 # Internal friction est
        spp = p_bit + total_loss + p_pipe

        # -----------------------------
        # ZERO FLOW THERMAL CALC
        # -----------------------------
        zf_mode = data.get('zero_flow_mode', 'immediate')
        zf_time = float(data.get('zero_flow_time', 12)) if zf_mode == 'insitu' else 0
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

        return jsonify({
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
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)