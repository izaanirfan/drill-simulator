import math

def minimum_curvature(p1, p2):
    md1, inc1, azi1 = p1
    md2, inc2, azi2 = p2

    inc1 = math.radians(inc1)
    inc2 = math.radians(inc2)
    azi1 = math.radians(azi1)
    azi2 = math.radians(azi2)

    dmd = md2 - md1

    dogleg = math.acos(
        math.cos(inc2 - inc1) -
        math.sin(inc1) * math.sin(inc2) *
        (1 - math.cos(azi2 - azi1))
    )

    rf = 1 if dogleg == 0 else 2 / dogleg * math.tan(dogleg / 2)

    dTVD = dmd / 2 * (math.cos(inc1) + math.cos(inc2)) * rf

    return dTVD

def build_trajectory(surveys):
    tvd = 0
    profile = []

    for i in range(1, len(surveys)):
        p1 = (surveys[i-1].md, surveys[i-1].inc, surveys[i-1].azi)
        p2 = (surveys[i].md, surveys[i].inc, surveys[i].azi)

        dTVD = minimum_curvature(p1, p2)
        tvd += dTVD

        profile.append({
            "md": surveys[i].md,
            "tvd": tvd,
            "inc": surveys[i].inc
        })

    return profile