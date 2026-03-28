import math

# -----------------------------
# AREA
# -----------------------------
def annular_area(outer_d, inner_d):
    return math.pi / 4 * ((outer_d / 12)**2 - (inner_d / 12)**2)

# -----------------------------
# HYD DIAMETER
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
# GET PIPE OD AT DEPTH
# -----------------------------
def get_pipe_od(bha_profile, md):

    for comp in bha_profile:
        if comp["top"] <= md <= comp["bottom"]:
            return comp["od"]

    return bha_profile[0]["od"]

# -----------------------------
# GET ANNULUS DIAMETER
# -----------------------------
def get_annulus_diameter(sections, md):

    for sec in sections:

        if sec.top_md <= md <= sec.end_md:

            t = sec.type.lower()

            if t == "open hole":
                return sec.hole_d

            if t == "casing":
                return sec.casing_id

            if t == "liner":

                if md >= sec.top_md:
                    return sec.casing_id

                # find parent casing
                for parent in sections:
                    if parent.type.lower() == "casing" and parent.top_md <= md <= parent.end_md:
                        return parent.casing_id

    # fallback
    last = sections[-1]
    return last.hole_d if last.hole_d > 0 else last.casing_id

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
# MAIN
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

        # ✅ NEW CORRECT LOGIC
        outer_d = get_annulus_diameter(data.well_sections, md2)
        pipe_od = get_pipe_od(bha_profile, md2)

        outer_d = max(outer_d, pipe_od + 0.1)

        area = annular_area(outer_d, pipe_od)
        area = max(area, 0.01)

        v = Q / (24.5 * area)

        dh = hydraulic_diameter(outer_d, pipe_od)
        dh = max(dh, 0.01)

        dp = pressure_loss(mw, v, dh, tau_y, K, n, length)
        cumulative_dp += dp

        ecd = mw + cumulative_dp / (0.051948 * md2)

        # Temperature correction
        temp = data.temperature
        surface_temp = temp.surface_temp
        bhct = temp.bhct

        temp_local = surface_temp + (bhct - surface_temp) * (md2 / depth)
        factor = 1 - 0.0003 * (temp_local - surface_temp)

        ecd_t = ecd * factor
        esd = mw * factor

        depths.append(md2)
        ecd_profile.append(round(ecd_t, 3))
        esd_profile.append(round(esd, 3))

    return {
        "summary": {
            "ecd_bottom": ecd_profile[-1],
            "esd_bottom": esd_profile[-1],
            "pressure_loss": round(cumulative_dp, 2)
        },
        "profile": {
            "depth": depths,
            "ecd": ecd_profile,
            "esd": esd_profile
        }
    }