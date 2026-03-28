import math

# -----------------------------
# CONSTANTS
# -----------------------------
GC = 32.174  # lbm-ft/lbf-s²

# -----------------------------
# AREA
# -----------------------------
def annular_area(Do, Di):
    return math.pi / 4 * ((Do/12)**2 - (Di/12)**2)

# -----------------------------
# HYD DIAMETER
# -----------------------------
def hydraulic_diameter(Do, Di):
    return (Do - Di) / 12


# -----------------------------
# BHA PROFILE
# -----------------------------
def build_bha_profile(bha, depth):

    profile = []
    current = depth

    for comp in reversed(bha):
        top = max(current - comp.length, 0)

        profile.append({
            "top": top,
            "bottom": current,
            "od": comp.od
        })

        current = top

    return profile


def get_pipe_od(profile, md):
    for p in profile:
        if p["top"] <= md <= p["bottom"]:
            return p["od"]
    return profile[-1]["od"]


# -----------------------------
# ANNULUS
# -----------------------------
def get_annulus(sections, md):

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

                for parent in sections:
                    if parent.type.lower() == "casing" and parent.top_md <= md <= parent.end_md:
                        return parent.casing_id

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

    n = math.log(tau1/tau2) / math.log(gamma1/gamma2)
    K = tau1 / (gamma1**n)

    return tau_y, K, n


# -----------------------------
# GENERALIZED REYNOLDS
# -----------------------------
def reynolds_hb(rho, v, dh, K, n):

    return (rho * (v**(2-n)) * (dh**n)) / (K * (8**(n-1)))


# -----------------------------
# FRICTION FACTOR
# -----------------------------
def friction_factor(Re):

    if Re < 2100:
        return 16 / Re

    # turbulent (approx)
    f = 0.005

    for _ in range(10):
        f = 1 / ( (4*math.log10(Re*math.sqrt(f)) - 0.4)**2 )

    return f


# -----------------------------
# PRESSURE LOSS
# -----------------------------
def pressure_loss(rho, v, dh, f, L):

    return f * (rho * v**2 / (2 * GC)) * (L / dh)


# -----------------------------
# MAIN
# -----------------------------
def run_simulation(data):

    Q = data.flowrate
    depth = data.depth
    mw = data.fluid.mw

    rho = mw * 8.34  # lbm/gal → lbm/ft³ approx

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

    step = 100
    md_list = list(range(step, int(depth)+step, step))

    for md in md_list:

        Do = get_annulus(data.well_sections, md)
        Di = get_pipe_od(bha_profile, md)

        Do = max(Do, Di + 0.1)

        area = annular_area(Do, Di)
        v = Q / (24.5 * area)

        dh = hydraulic_diameter(Do, Di)

        Re = reynolds_hb(rho, v, dh, K, n)
        f = friction_factor(Re)

        dp = pressure_loss(rho, v, dh, f, step)
        cumulative_dp += dp

        ecd = mw + cumulative_dp / (0.051948 * md)

        # temperature
        Tsurf = data.temperature.surface_temp
        Tbh = data.temperature.bhct

        T = Tsurf + (Tbh - Tsurf)*(md/depth)
        factor = 1 - 0.0003*(T - Tsurf)

        ecd_t = ecd * factor
        esd = mw * factor

        depths.append(md)
        ecd_profile.append(round(ecd_t,3))
        esd_profile.append(round(esd,3))

    return {
        "summary": {
            "ecd_bottom": ecd_profile[-1],
            "esd_bottom": esd_profile[-1]
        },
        "profile": {
            "depth": depths,
            "ecd": ecd_profile,
            "esd": esd_profile
        }
    }