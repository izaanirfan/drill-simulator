import math

# -----------------------------
# CONSTANTS
# -----------------------------
GC = 32.174  # gravitational constant (lbm-ft/lbf-s²)

# -----------------------------
# GEOMETRY
# -----------------------------
def annular_area(Do, Di):
    # inches → ft
    Do_ft = Do / 12
    Di_ft = Di / 12
    area = math.pi / 4 * (Do_ft**2 - Di_ft**2)
    return max(area, 0.0001)


def hydraulic_diameter(Do, Di):
    return max((Do - Di) / 12, 0.01)  # ft


# -----------------------------
# BHA PROFILE
# -----------------------------
def build_bha_profile(bha, depth):

    profile = []
    current = depth

    if not bha:
        return [{"top": 0, "bottom": depth, "od": 5}]

    for comp in reversed(bha):
        length = comp.length or 0
        od = comp.od or 5

        top = max(current - length, 0)

        profile.append({
            "top": top,
            "bottom": current,
            "od": od
        })

        current = top

    return profile


def get_pipe_od(profile, md):
    for p in profile:
        if p["top"] <= md <= p["bottom"]:
            return p["od"]
    return profile[-1]["od"]


# -----------------------------
# WELL GEOMETRY
# -----------------------------
def get_annulus(sections, md):

    if not sections:
        return 8.5

    for sec in sections:

        if sec.top_md <= md <= sec.end_md:

            t = sec.type.lower()

            if t == "open hole":
                return sec.hole_d or 8.5

            if t == "casing":
                return sec.casing_id or 8.5

            if t == "liner":

                if md >= sec.top_md:
                    return sec.casing_id or 8.5

                # parent casing
                for parent in sections:
                    if parent.type.lower() == "casing" and parent.top_md <= md <= parent.end_md:
                        return parent.casing_id or 8.5

    return sections[-1].hole_d or 8.5


# -----------------------------
# HB FIT
# -----------------------------
def fit_hb(f600, f300, f3):

    tau_y = f3  # lb/100 ft²

    gamma1 = 1022
    gamma2 = 511

    tau1 = max(f600 - tau_y, 0.1)
    tau2 = max(f300 - tau_y, 0.1)

    n = math.log(tau1 / tau2) / math.log(gamma1 / gamma2)
    K = tau1 / (gamma1 ** n)

    return tau_y, K, n


# -----------------------------
# HB → APPARENT VISCOSITY
# -----------------------------
def apparent_viscosity_hb(tau_y, K, n, v, dh):

    # shear rate
    gamma = max(8 * v / dh, 0.1)

    # shear stress
    tau = tau_y + K * (gamma ** n)

    # viscosity (lb·s/ft²)
    mu = tau / gamma

    # convert to cp
    mu_cp = mu * 4788

    return max(mu_cp, 1)


# -----------------------------
# MAIN SIMULATION
# -----------------------------
def run_simulation(data):

    try:

        Q = float(data.flowrate or 400)      # gpm
        depth = max(float(data.depth or 10000), 100)
        mw = float(data.fluid.mw or 10)

        # density (lbm/ft³)
        rho = mw * 0.052 * 144

        # HB rheology
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

        step = max(int(depth / 100), 50)

        md_list = list(range(step, int(depth) + step, step))

        for md in md_list:

            # geometry
            Do = get_annulus(data.well_sections, md)
            Di = get_pipe_od(bha_profile, md)

            Do = max(Do, Di + 0.1)

            area = annular_area(Do, Di)

            # velocity ft/s
            v = Q / (24.5 * area)

            dh = hydraulic_diameter(Do, Di)

            # HB apparent viscosity
            mu_cp = apparent_viscosity_hb(tau_y, K, n, v, dh)

            # Reynolds number (field form)
            Re = (928 * rho * v * dh) / mu_cp

            # friction factor
            if Re < 2100:
                f = 16 / max(Re, 1)
            else:
                f = 0.079 / (Re ** 0.25)

            # pressure loss (psi)
            dp = f * (rho * v**2 / (2 * GC)) * (step / dh)

            cumulative_dp += dp

            # ECD
            ecd = mw + cumulative_dp / (0.051948 * md)

            # temperature correction
            Tsurf = data.temperature.surface_temp
            Tbh = data.temperature.bhct

            T = Tsurf + (Tbh - Tsurf) * (md / depth)
            factor = 1 - 0.0003 * (T - Tsurf)

            ecd_t = ecd * factor
            esd = mw * factor

            depths.append(md)
            ecd_profile.append(round(ecd_t, 3))
            esd_profile.append(round(esd, 3))

        return {
            "summary": {
                "ecd_bottom": ecd_profile[-1] if ecd_profile else 0,
                "esd_bottom": esd_profile[-1] if esd_profile else 0,
                "pressure_loss": round(cumulative_dp, 2)
            },
            "profile": {
                "depth": depths,
                "ecd": ecd_profile,
                "esd": esd_profile
            }
        }

    except Exception as e:
        return {"error": str(e)}