import math
from rheology import hb_fit
from temperature import *
from trajectory import build_trajectory
from cuttings import *

GPM_TO_FT3S = 0.002228
IN2_TO_FT2 = 1 / 144

def annular_area(hole_d, pipe_od):
    return math.pi / 4 * (hole_d**2 - pipe_od**2) * IN2_TO_FT2

def annular_velocity(q_gpm, hole_d, pipe_od):
    q = q_gpm * GPM_TO_FT3S
    area = annular_area(hole_d, pipe_od)
    return q / area

def hydraulic_diameter(hole_d, pipe_od):
    return (hole_d - pipe_od) / 12

def reynolds_hb(rho, v, Dh, K, n):
    return (rho * v**(2-n) * Dh**n) / (K * 8**(n-1))

def friction_factor(Re):
    return 16/Re if Re < 2100 else 0.079/(Re**0.25)

def segment_pressure_loss(f, L, Dh, rho, v):
    return f * (L/Dh) * (rho * v**2 / 2) / 144

def run_simulation(data):

    fann = [
        data.fluid.fann_600,
        data.fluid.fann_300,
        data.fluid.fann_200,
        data.fluid.fann_100,
        data.fluid.fann_6,
        data.fluid.fann_3
    ]

    tau0, K, n = hb_fit(fann)

    traj = build_trajectory(data.trajectory)

    ecd_profile = []
    depth_profile = []
    transport_ratio_profile = []
    bed_height_profile = []

    cumulative_annular = 0

    for point in traj:

        md = point["md"]
        tvd = point["tvd"]
        inc = point["inc"]

        T_trans = temp_transient(tvd, data.temperature.surface_temp, data.temperature.bhct, data.depth)
        mw = density_corrected(data.fluid.mw, data.temperature.beta, T_trans)

        rho = mw * 8.34

        for comp in data.bha:
            v = annular_velocity(data.flowrate, data.hole_diameter, comp.od)
            Dh = hydraulic_diameter(data.hole_diameter, comp.od)

            Re = reynolds_hb(rho, v, Dh, K, n)
            f = friction_factor(Re)

            dp = segment_pressure_loss(f, comp.length, Dh, rho, v)

            correction = 1 + 0.4 * math.sin(math.radians(inc))
            cumulative_annular += dp * correction

        hydro = 0.052 * mw * tvd
        bhp = hydro + cumulative_annular + data.sbp
        ecd = bhp / (0.052 * tvd)

        # Cuttings
        rho_cut = data.cuttings.density * 8.34
        Vs = slip_velocity(data.cuttings.size, rho_cut, rho)
        TR = v / max(Vs, 0.01)
        TR_eff = TR * (1 - 0.5 * math.sin(math.radians(inc)))
        bed = bed_height(TR_eff, inc)

        depth_profile.append(tvd)
        ecd_profile.append(ecd)
        transport_ratio_profile.append(TR_eff)
        bed_height_profile.append(bed)

    return {
        "depth_profile": depth_profile,
        "ecd_profile": ecd_profile,
        "transport_ratio": transport_ratio_profile,
        "bed_height": bed_height_profile
    }