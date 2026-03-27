import math

def slip_velocity(size, rho_c, rho_m):
    g = 32.2
    d = size / 12
    return 1.5 * ((g * d * (rho_c - rho_m) / rho_m) ** 0.5)

def cuttings_concentration(rop, hole_d, flowrate):
    area = math.pi/4 * hole_d**2
    q = flowrate * 0.002228
    return (rop / 60) * area / q

def bed_height(TR_eff, inc):
    if inc < 30:
        return 0
    k = 0.3
    return k / max(TR_eff, 0.1)