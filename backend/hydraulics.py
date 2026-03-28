import math

# -----------------------------
# GEOMETRY
# -----------------------------
def annular_area(Do_in, Di_in):
    """Returns area in square inches."""
    return max((math.pi / 4) * (Do_in**2 - Di_in**2), 0.1)

def hydraulic_diameter(Do_in, Di_in):
    """Returns hydraulic diameter (dh) in inches."""
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
    return 5.0  # Fallback to standard Drill Pipe OD

def get_annulus_id(sections, md):
    for sec in sections:
        if sec.top_md <= md <= sec.end_md:
            t = sec.type.lower()
            return sec.hole_d if t == "open hole" else sec.casing_id
    return sections[-1].hole_d

# -----------------------------
# RHEOLOGY (Herschel-Bulkley)
# -----------------------------
def fit_hb(f600, f300, f3):
    tau_y = f3 
    # n = log(tau1/tau2) / log(gamma1/gamma2)
    n = math.log10((f600 - tau_y) / (f300 - tau_y)) / math.log10(600 / 300)
    K = (f300 - tau_y) / (511**n)
    return tau_y, K, n

# -----------------------------
# REFINED PRESSURE LOSS MODEL
# -----------------------------
def calculate_annular_pressure_loss(mw, flowrate, Do, Di, mu_cp, length):
    """
    Calculates annular pressure loss in psi.
    mw: ppg, flowrate: gpm, Do/Di: inches, length: feet
    """
    dh = Do - Di
    # Annular Velocity in ft/min
    v = (24.48 * flowrate) / (Do**2 - Di**2)
    
    # Reynolds Number (Generalised for non-Newtonian fluids)
    reynolds = (15.47 * mw * v * dh) / mu_cp
    
    # Friction Factor (f)
    if reynolds < 2100:
        f = 64 / reynolds # Laminar
    else:
        f = 0.3164 / (reynolds**0.25) # Turbulent (Blasius)
    
    # Pressure loss (psi) using dP = (f * L * rho * v^2) / (25.6 * dh)
    # Note: v converted from ft/min to ft/sec for the standard Darcy-Weisbach form
    v_sec = v / 60
    dp = (f * length * mw * (v_sec**2)) / (25.6 * (dh / 12)) 
    
    # Convert result to psi based on standard field units
    # Simplified field version of the above:
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
        precision_const = 0.051948

        tau_y, K, n = fit_hb(data.fluid.fann_600, data.fluid.fann_300, data.fluid.fann_3)
        bha_profile = build_bha_profile(data.bha, total_depth)

        depths, ecd_profile, esd_profile = [], [], []
        cumulative_annular_loss = 0
        
        # Determine calculation resolution
        step = max(int(total_depth / 100), 50)
        md_points = list(range(step, int(total_depth) + step, step))

        for md in md_points:
            if md > total_depth: md = total_depth

            Do = get_annulus_id(data.well_sections, md)
            Di = get_pipe_od(bha_profile, md)
            dh = hydraulic_diameter(Do, Di)
            
            # 1. Calculate Shear Rate (gamma) and Apparent Viscosity (mu_cp)
            v_ftmin = (24.48 * Q) / (Do**2 - Di**2)
            gamma = (1.6 * v_ftmin) / dh 
            mu_cp = ((tau_y + K * (gamma**n)) / gamma) * 478.8 if gamma > 0.1 else 100
            
            # 2. Calculate dP for this specific segment length
            dp = calculate_annular_pressure_loss(mw, Q, Do, Di, mu_cp, step)
            cumulative_annular_loss += dp

            # 3. ECD = MW + (Sum of dP / (Constant * TVD))
            ecd_base = mw + (cumulative_annular_loss / (precision_const * md))

            # 4. Temperature Correction (ESD and ECD)
            Tsurf = data.temperature.surface_temp
            Tbh = data.temperature.bhct
            T_current = Tsurf + (Tbh - Tsurf) * (md / total_depth)
            temp_factor = 1 - 0.0003 * (T_current - Tsurf)

            depths.append(md)
            ecd_profile.append(round(ecd_base * temp_factor, 3))
            esd_profile.append(round(mw * temp_factor, 3))

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