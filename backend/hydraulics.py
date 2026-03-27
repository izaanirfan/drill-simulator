import math

def run_simulation(data):

    # -----------------------------
    # BASIC INPUTS
    # -----------------------------
    Q = data.flowrate            # gpm
    depth = data.depth           # ft
    hole_d = data.hole_diameter  # inch

    # -----------------------------
    # FLUID PROPERTIES
    # -----------------------------
    mw = data.fluid.mw  # ppg

    # Convert MW to density (ppg → lb/ft3)
    rho_fluid = mw * 7.48

    # Fann readings
    f600 = data.fluid.fann_600
    f300 = data.fluid.fann_300

    # Simple PV/YP estimation (Bingham approx)
    pv = f600 - f300
    yp = f300 - pv

    # -----------------------------
    # GEOMETRY
    # -----------------------------
    hole_area = math.pi * (hole_d / 12)**2 / 4   # ft2

    # Velocity (ft/s)
    v = Q / (24.5 * hole_area)

    # -----------------------------
    # PRESSURE LOSS (VERY SIMPLIFIED)
    # -----------------------------
    friction_factor = 0.02  # placeholder

    dp = friction_factor * depth * mw

    # -----------------------------
    # ECD CALCULATION
    # -----------------------------
    ecd = mw + (dp / (0.052 * depth))

    # -----------------------------
    # CUTTINGS TRANSPORT
    # -----------------------------
    cut = data.cuttings

    rho_cuttings = cut.density * 62.4   # convert SG → lb/ft3
    d_cut = cut.size / 12              # inch → ft
    rpm = cut.rpm

    g = 32.2

    # SAFE SETTLING VELOCITY
    value = (rho_cuttings - rho_fluid) * g * d_cut

    # Prevent negative → no complex numbers
    value = max(value, 0)

    Vs = math.sqrt(value)

    # Ensure real and safe
    if isinstance(Vs, complex):
        Vs = abs(Vs)

    Vs = max(Vs, 0.01)

    # Transport Ratio
    TR = v / Vs

    # Clamp TR for realism
    TR = min(TR, 5)

    # -----------------------------
    # TEMPERATURE EFFECT (SIMPLIFIED)
    # -----------------------------
    temp = data.temperature

    surface_temp = temp.surface_temp
    bhct = temp.bhct

    temp_factor = 1 - 0.0003 * (bhct - surface_temp)

    ecd_temp_corrected = ecd * temp_factor

    # -----------------------------
    # OUTPUT PROFILE (SIMPLE)
    # -----------------------------
    depths = []
    ecd_profile = []

    step = depth / 20

    for i in range(21):
        d = i * step
        depths.append(d)
        ecd_profile.append(ecd_temp_corrected)

    # -----------------------------
    # RESULTS
    # -----------------------------
    return {
        "ecd_surface": round(ecd, 2),
        "ecd_temp_corrected": round(ecd_temp_corrected, 2),
        "pressure_loss": round(dp, 2),
        "transport_ratio": round(TR, 2),
        "velocity": round(v, 2),
        "settling_velocity": round(Vs, 2),
        "profile": {
            "depth": depths,
            "ecd": ecd_profile
        }
    }