import math

# -----------------------------
# ANNULAR AREA
# -----------------------------
def annular_area(outer_d, inner_d):
    return math.pi / 4 * ((outer_d / 12)**2 - (inner_d / 12)**2)

# -----------------------------
# HYDRAULIC DIAMETER
# -----------------------------
def hydraulic_diameter(outer_d, inner_d):
    return (outer_d - inner_d) / 12

# -----------------------------
# BUILD BHA PROFILE
# -----------------------------
def build_bha_profile(bha, total_depth):

    profile = []
    current_depth = total_depth

    for comp in reversed(bha):
        top = current_depth - comp.length

        profile.append({
            "top": top,
            "bottom": current_depth,
            "od": comp.od
        })

        current_depth = top

    return profile

# -----------------------------
# GET PIPE OD
# -----------------------------
def get_pipe_od(bha_profile, md):
    for comp in bha_profile:
        if comp["top"] <= md <= comp["bottom"]:
            return comp["od"]
    return bha_profile[0]["od"]

# -----------------------------
# GET ANNULUS DIAMETER (NEW LOGIC)
# -----------------------------
def get_annulus_diameter(sections, md):

    current_section = None

    for sec in sections:
        if sec.top_md <= md <= sec.end_md:
            current_section = sec
            break

    if current_section is None:
        return sections[-1].hole_d or sections[-1].casing_id

    # OPEN HOLE
    if current_section.type.lower() == "open hole":
        return current_section.hole_d

    # CASING
    if current_section.type.lower() == "casing":
        return current_section.casing_id

    # LINER
    if current_section.type.lower() == "liner":

        if md >= current_section.liner_top:
            return current_section.casing_id
        else:
            # find parent casing
            for sec in sections:
                if sec.type.lower() == "casing" and sec.top_md <= md <= sec.end_md:
                    return sec.casing_id

    # fallback
    return current_section.casing_id or current_section.hole_d

# -----------------------------
# HB FIT
# -----------------------------
def fit_hb(f600, f300, f3):

    tau_y = f3

    gamma1 = 1022
    gamma2 = 511

    tau1 = max(f600 - tau_y, 0.1)
    tau2 = max(f300 - tau_y, 0.1)

    n = math.log(tau1 / tau2) / math.log(gamma1 / gamma2)
    K = tau1 / (gamma1 ** n)

    return tau_y, K, n

# -----------------------------
# APPARENT VISCOSITY
# -----------------------------
def apparent_viscosity(tau_y, K, n, gamma):
    tau = tau_y + K * (gamma ** n)
    return tau / max(gamma, 0.1)

# -----------------------------
# PRESSURE LOSS
# -----------------------------
def pressure_loss(mw, v, dh, tau_y, K, n, length):

    gamma = 8 * v / dh
    mu = apparent_viscosity(tau_y, K, n, gamma)

    rho = mw * 7.48
    Re = (rho * v * dh) / max(mu, 0.01)

    if Re < 2100:
        dp = (32 * mu * v / (dh**2)) * length * 0.01
    else:
        f = 0.079 / (Re ** 0.25)
        dp = f * (rho * v**2 / (2 * dh)) * length * 0.01

    return dp

# -----------------------------
# MAIN SIMULATION
# -----------------------------
def run_simulation(data):

    Q = data.flowrate
    depth = data.depth
    mw = data.fluid.mw

    tau_y, K, n = fit_hb(
        data.fluid.fann_600,
        data.fluid.fann_300,
        data.fluid.fann_3
    )

    bha_profile = build_bha_profile(data.bha, depth)

    depths = []
    ecd_profile = []
    esd_profile = []

    cumulative_dp = 0

    for i in range(1, len(data.trajectory)):

        md1 = data.trajectory[i-1].md
        md2 = data.trajectory[i].md
        length = md2 - md1

        # 🔥 NEW CORE LOGIC
        outer_d = get_annulus_diameter(data.well_sections, md2)
        pipe_od = get_pipe_od(bha_profile, md2)

        area = annular_area(outer_d, pipe_od)
        area = max(area, 0.01)

        v = Q / (24.5 * area)

        dh = hydraulic_diameter(outer_d, pipe_od)
        dh = max(dh, 0.01)

        dp = pressure_loss(mw, v, dh, tau_y, K, n, length)

        cumulative_dp += dp

        ecd = mw + cumulative_dp / (0.051948 * md2)

        temp = data.temperature
        temp_local = temp.surface_temp + (temp.bhct - temp.surface_temp) * (md2 / depth)
        temp_factor = 1 - 0.0003 * (temp_local - temp.surface_temp)

        ecd_temp = ecd * temp_factor
        esd = mw * temp_factor

        depths.append(md2)
        ecd_profile.append(round(ecd_temp, 3))
        esd_profile.append(round(esd, 3))

    return {
        "summary": {
            "ecd_bottom": ecd_profile[-1],
            "esd_bottom": esd_profile[-1],
            "tau_y": round(tau_y, 2),
            "K": round(K, 4),
            "n": round(n, 3)
        },
        "profile": {
            "depth": depths,
            "ecd": ecd_profile,
            "esd": esd_profile
        }
    }