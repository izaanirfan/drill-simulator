import math

# -----------------------------
# ANNULAR AREA
# -----------------------------
def annular_area(hole_d, pipe_od):
    return math.pi / 4 * ((hole_d / 12)**2 - (pipe_od / 12)**2)

# -----------------------------
# HYDRAULIC DIAMETER
# -----------------------------
def hydraulic_diameter(hole_d, pipe_od):
    return (hole_d - pipe_od) / 12

# -----------------------------
# BUILD BHA PROFILE (TOP → BOTTOM)
# -----------------------------
def build_bha_profile(bha, total_depth):

    profile = []
    current_depth = total_depth

    # Build from bottom up
    for comp in reversed(bha):
        length = comp.length
        top = current_depth - length

        profile.append({
            "top": top,
            "bottom": current_depth,
            "od": comp.od
        })

        current_depth = top

    return profile

# -----------------------------
# GET PIPE OD AT DEPTH
# -----------------------------
def get_pipe_od(bha_profile, md):

    for comp in bha_profile:
        if comp["top"] <= md <= comp["bottom"]:
            return comp["od"]

    return bha_profile[0]["od"]  # fallback

# -----------------------------
# GET HOLE SIZE
# -----------------------------
def get_hole_size(sections, md):

    for sec in sections:
        if sec.top_md <= md <= sec.end_md:
            return sec.hole_d

    return sections[-1].hole_d

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

    f600 = data.fluid.fann_600
    f300 = data.fluid.fann_300
    f3 = data.fluid.fann_3

    tau_y, K, n = fit_hb(f600, f300, f3)

    sections = data.well_sections
    bha = data.bha
    traj = data.trajectory

    bha_profile = build_bha_profile(bha, depth)

    depths = []
    ecd_profile = []
    esd_profile = []

    cumulative_dp = 0

    for i in range(1, len(traj)):

        md1 = traj[i-1].md
        md2 = traj[i].md
        length = md2 - md1

        # 🔥 NEW LOGIC
        hole_d = get_hole_size(sections, md2)
        pipe_od = get_pipe_od(bha_profile, md2)

        area = annular_area(hole_d, pipe_od)
        area = max(area, 0.01)

        v = Q / (24.5 * area)

        dh = hydraulic_diameter(hole_d, pipe_od)
        dh = max(dh, 0.01)

        dp = pressure_loss(mw, v, dh, tau_y, K, n, length)

        cumulative_dp += dp

        ecd = mw + cumulative_dp / (0.051948 * md2)

        temp = data.temperature
        surface_temp = temp.surface_temp
        bhct = temp.bhct

        temp_local = surface_temp + (bhct - surface_temp) * (md2 / depth)
        temp_factor = 1 - 0.0003 * (temp_local - surface_temp)

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