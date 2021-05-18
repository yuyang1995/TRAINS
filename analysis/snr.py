import numpy as np
from param import parset
from residual import FullResiduals

if __name__ == "__main__":
    alphaS, sin_deltaS, cos_iota, psi, phi0, log_zeta, log_omega, log_tm \
        = np.loadtxt("../data/sim_source.txt", unpack=True, usecols=range(8))
    zeta, omega, tm = 10**log_zeta, 10**log_omega, 10**log_tm

    alphaP, deltaP, distP, sigma = np.loadtxt(
        "../{:s}".format(parset['PulsarCatalog']),
        unpack=True, usecols=[0, 1, 2, 4])
    Np = alphaP.size

    t = np.loadtxt("../data/sim_ptr.txt", usecols=0)
    Nt = t.size // Np
    t = t.reshape((Np, Nt))

    flag_evolve = int(parset['FlagEvolve'])
    res = FullResiduals(alphaS, sin_deltaS, cos_iota, psi, phi0,
                        zeta, omega, tm, alphaP, deltaP, t, distP,
                        None, 1, flag_evolve)

    rhoI = np.linalg.norm(res, axis=-1) / sigma[:, None]
    rho = np.linalg.norm(rhoI, axis=0)
    np.savetxt("../data/sim_snr.txt", rho, fmt='%.6e')

    signal = np.einsum("ijk->ik", res, optimize=True)
    np.savetxt("../data/sim_signal.txt", signal.flatten(), fmt='%.6e')

    # noise = np.random.normal(size=signal.size).reshape(signal.shape) * 15e-7
    # data = np.stack((t.flatten(), (signal + noise).flatten()), axis=-1)
    # np.savetxt("../data/sim_ptr.txt", data, fmt='%.6e')

    print(rho.argsort()[-10:])
    print(np.sort(rho)[-10:])
    print(log_omega[rho.argsort()][-10:])
    if flag_evolve:
        print(log_tm[rho.argsort()][-10:])
