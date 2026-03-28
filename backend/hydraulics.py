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
# FIT HB PARAMETERS
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
# REYNOLDS NUMBER
# -----------------------------
def reynolds_number(mw, v, dh, mu):
    rho = mw * 7.48
    return (rho * v * dh) / max(mu, 0.01)

# -----------------------------
# PRESSURE LOSS
# -----------------------------
def pressure_loss(mw, v, dh, tau_y, K, n, length):

    gamma = 8 * v / dh
    mu = apparent_viscosity(tau_y, K, n, gamma)
    Re = reynolds_number(mw, v, dh, mu)

    if Re < 2100:
        dp = (32 * mu * v / (dh**2)) * length * 0.01
    else:
        f = 0.079 / (Re ** 0.25)
        rho = mw * 7.48
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
    traj = data.trajectory

    depths = []
    ecd_profile = []
    esd_profile = []

    cumulative_dp = 0

    for i in range(1, len(traj)):

        md1 = traj[i-1].md
        md2 = traj[i].md
        length = md2 - md1

        # -----------------------------
        # FIND SECTION (FIXED)
        # -----------------------------
        hole_d = None
        pipe_od = None

        for sec in sections:
            if md2 <= sec.end_md:   # ✅ FIX HERE
                hole_d = sec.hole_d
                pipe_od = sec.pipe_od
                break

        if hole_d is None:
            hole_d = sections[-1].hole_d
            pipe_od = sections[-1].pipe_od

        # -----------------------------
        # CALCULATIONS
        # -----------------------------
        area = annular_area(hole_d, pipe_od)
        area = max(area, 0.01)

        v = Q / (24.5 * area)

        dh = hydraulic_diameter(hole_d, pipe_od)
        dh = max(dh, 0.01)

        dp = pressure_loss(mw, v, dh, tau_y, K, n, length)
        cumulative_dp += dp

        # ECD (correct constant)
        ecd = mw + cumulative_dp / (0.051948 * md2)

        # Temperature
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
            "n": round(n, 3),
            "total_pressure_loss": round(cumulative_dp, 2)
        },
        "profile": {
            "depth": depths,
            "ecd": ecd_profile,
            "esd": esd_profile
        }
    }