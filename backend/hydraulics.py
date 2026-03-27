def run_simulation(data):

    flowrate = data.flowrate
    depth = data.depth
    mw = data.fluid.mw

    # -----------------------------
    # SIMPLE HYDRAULICS MODEL
    # -----------------------------

    pressure_loss = 0.02 * flowrate * depth / 100  # simplified

    # Temperature correction (simple)
    temp_factor = 0.98

    # ECD
    ecd_bottom = mw + pressure_loss / (0.052 * depth)
    ecd_temp = ecd_bottom * temp_factor

    # ESD (no friction)
    esd_bottom = mw * temp_factor

    # -----------------------------
    # PROFILE (DEPTH VARIATION)
    # -----------------------------
    depth_list = []
    ecd_list = []
    esd_list = []

    steps = 20

    for i in range(steps + 1):
        d = depth * i / steps

        # simulate higher losses at bottom (smaller annulus)
        factor = 1 + 0.3 * (d / depth)

        ecd_d = mw + (pressure_loss * factor) / (0.052 * depth)
        ecd_d = ecd_d * temp_factor

        esd_d = mw * temp_factor

        depth_list.append(d)
        ecd_list.append(ecd_d)
        esd_list.append(esd_d)

    # -----------------------------
    # FINAL OUTPUT (FIXED FORMAT)
    # -----------------------------
    return {
        "summary": {
            "ecd_bottom": round(ecd_temp, 2),
            "esd_bottom": round(esd_bottom, 2)
        },
        "profile": {
            "depth": depth_list,
            "ecd": ecd_list,
            "esd": esd_list
        }
    }