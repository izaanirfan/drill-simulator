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

    # 🔥 KEY FIX: create depth steps
    step = 200   # ft (you can change later)

    md_list = list(range(step, int(depth)+step, step))

    for md2 in md_list:

        length = step

        # ANNULUS
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
            "esd_bottom": esd_profile[-1]
        },
        "profile": {
            "depth": depths,
            "ecd": ecd_profile,
            "esd": esd_profile
        }
    }