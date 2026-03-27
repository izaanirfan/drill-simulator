import numpy as np

def hb_fit(fann):
    rpm = np.array([600,300,200,100,6,3])
    shear_rate = 1.703 * rpm
    shear_stress = np.array(fann)

    log_gamma = np.log(shear_rate[1:])
    log_tau = np.log(shear_stress[1:])

    n, logK = np.polyfit(log_gamma, log_tau, 1)
    K = np.exp(logK)

    tau0 = max(0, shear_stress[-1] - K*(shear_rate[-1]**n))

    return tau0, K, n