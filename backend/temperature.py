def temp_transient(depth, T_in, bhct, total_depth):
    return T_in + (bhct - T_in) * (depth / total_depth)

def temp_geothermal(depth, surface_temp, grad):
    return surface_temp + grad * depth

def density_corrected(mw, beta, T, T_ref=80):
    return mw * (1 - beta * (T - T_ref))

def rheology_temp_correction(K, tau0, T, T_ref=80):
    import math
    a = 0.02
    b = 0.015

    K_corr = K * math.exp(-a * (T - T_ref))
    tau0_corr = tau0 * math.exp(-b * (T - T_ref))

    return K_corr, tau0_corr