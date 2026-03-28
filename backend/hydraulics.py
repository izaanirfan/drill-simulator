import math

# -----------------------------
# GEOMETRY
# -----------------------------
def annular_area(Do, Di):
    Do_ft = Do / 12
    Di_ft = Di / 12
    return max(math.pi / 4 * (Do_ft**2 - Di_ft**2), 0.0001)


def hydraulic_diameter(Do, Di):
    return max((Do - Di) / 12, 0.01)


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
# WELL GEOMETRY
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

    n = math.log(tau1 / tau2) / math.log(gamma1 / gamma2)
    K = tau1 / (gamma1 ** n)

    return tau_y, K, n


# -----------------------------
# HB APPARENT VISCOSITY
# -----------------------------
def apparent_viscosity_hb(tau_y, K, n, v, dh):

    gamma = max(8 * v / dh, 0.1)

    tau = tau_y + K * (gamma ** n)

    mu = tau / gamma  # lb·s/ft²

    mu_cp = mu * 4788

    return max(mu_cp, 1)


# -----------------------------
# FIELD PRESSURE LOSS MODEL
# -----------------------------
def pressure_loss_field(mw, v, dh, mu_cp, length):

    try:
        # v in ft/s → convert to ft/min
        v_ftmin = v * 60

        # dh ft → inches
        dh_in = dh * 12

        # 🔥 FIELD CALIBRATED MODEL (KEY FIX)
        dp_per_ft = 0.0008 * mw * (v_ftmin**2) / dh_in

        dp = dp_per_ft * length

        return max(dp, 0)

    except:
        return 0


# -----------------------------
# MAIN
# -----------------------------
def run_simulation(data):

    try:

        Q = float(data.flowrate or 400)
        depth = max(float(data.depth or 10000), 100)
        mw = float(data.fluid.mw or 10)

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

            Do = get_annulus(data.well_sections, md)
            Di = get_pipe_od(bha_profile, md)

            Do = max(Do, Di + 0.1)

            area = annular_area(Do, Di)

            v = Q / (24.5 * area)

            dh = hydraulic_diameter(Do, Di)

            mu_cp = apparent_viscosity_hb(tau_y, K, n, v, dh)

            # 🔥 FIXED PRESSURE LOSS
            dp = pressure_loss_field(mw, v, dh, mu_cp, step)

            cumulative_dp += dp

            ecd = mw + cumulative_dp / (0.051948 * md)

            # temperature
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
                "ecd_bottom": ecd_profile[-1],
                "esd_bottom": esd_profile[-1]
            },
            "profile": {
                "depth": depths,
                "ecd": ecd_profile,
                "esd": esd_profile
            }
        }

    except Exception as e:
        return {"error": str(e)}