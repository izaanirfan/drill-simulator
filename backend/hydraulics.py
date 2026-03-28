import math

# -----------------------------
# SAFE AREA
# -----------------------------
def annular_area(outer_d, inner_d):
    try:
        return max(math.pi / 4 * ((outer_d/12)**2 - (inner_d/12)**2), 0.01)
    except:
        return 0.01


def hydraulic_diameter(outer_d, inner_d):
    try:
        return max((outer_d - inner_d) / 12, 0.01)
    except:
        return 0.01


# -----------------------------
# BHA PROFILE
# -----------------------------
def build_bha_profile(bha, total_depth):

    profile = []
    current_depth = total_depth

    if not bha:
        return [{"top": 0, "bottom": total_depth, "od": 5}]

    for comp in reversed(bha):

        length = comp.length if comp.length else 0
        od = comp.od if comp.od else 5

        top = max(current_depth - length, 0)

        profile.append({
            "top": top,
            "bottom": current_depth,
            "od": od
        })

        current_depth = top

    return profile


def get_pipe_od(bha_profile, md):

    for comp in bha_profile:
        if comp["top"] <= md <= comp["bottom"]:
            return comp["od"]

    return bha_profile[-1]["od"] if bha_profile else 5


# -----------------------------
# ANNULUS
# -----------------------------
def get_annulus_diameter(sections, md):

    if not sections:
        return 8.5

    for sec in sections:

        try:
            if sec.top_md <= md <= sec.end_md:

                t = sec.type.lower()

                if t == "open hole":
                    return sec.hole_d or 8.5

                if t == "casing":
                    return sec.casing_id or 8.5

                if t == "liner":

                    if md >= sec.top_md:
                        return sec.casing_id or 8.5

                    for parent in sections:
                        if parent.type.lower() == "casing" and parent.top_md <= md <= parent.end_md:
                            return parent.casing_id or 8.5
        except:
            continue

    last = sections[-1]
    return last.hole_d if last.hole_d else 8.5


# -----------------------------
# HB FIT
# -----------------------------
def fit_hb(f600, f300, f3):

    try:
        tau_y = f3

        gamma1 = 1022
        gamma2 = 511

        tau1 = max(f600 - tau_y, 0.1)
        tau2 = max(f300 - tau_y, 0.1)

        n = math.log(tau1 / tau2) / math.log(gamma1 / gamma2)
        K = tau1 / (gamma1 ** n)

        return tau_y, K, n
    except:
        return 5, 0.5, 0.8


def apparent_viscosity(tau_y, K, n, gamma):

    gamma = max(gamma, 0.1)

    tau = tau_y + K * (gamma ** n)

    return tau / gamma


def pressure_loss(mw, v, dh, tau_y, K, n, length):

    try:
        # Apparent viscosity (cp)
        gamma = 8 * v / dh
        mu = apparent_viscosity(tau_y, K, n, gamma)

        # 🔥 FIELD CALIBRATED MODEL
        # dp (psi) ≈ (0.0001 * MW * v^2 / dh) * length
        dp = 0.0001 * mw * (v**2) / dh * length

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

        for md2 in md_list:

            if md2 <= 0:
                continue

            length = step

            outer_d = get_annulus_diameter(data.well_sections, md2)
            pipe_od = get_pipe_od(bha_profile, md2)

            outer_d = max(outer_d, pipe_od + 0.1)

            area = annular_area(outer_d, pipe_od)

            v = Q / (24.5 * area)

            dh = hydraulic_diameter(outer_d, pipe_od)

            dp = pressure_loss(mw, v, dh, tau_y, K, n, length)
            cumulative_dp += dp

            ecd = mw + cumulative_dp / (0.051948 * md2)

            # Temperature
            surface_temp = data.temperature.surface_temp
            bhct = data.temperature.bhct

            temp_local = surface_temp + (bhct - surface_temp) * (md2 / depth)
            factor = 1 - 0.0003 * (temp_local - surface_temp)

            ecd_t = ecd * factor
            esd = mw * factor

            depths.append(md2)
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
        return {
            "error": str(e)
        }