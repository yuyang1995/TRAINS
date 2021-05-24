import numpy as np
from param import parset
from plot_toolbox import resplot, skymap

if __name__ == "__main__":
    pulsars = np.loadtxt("../{:s}".format(parset['PulsarCatalog']))
    sources = np.loadtxt("../data/sim_source.txt", usecols=[0, 1])
    snr = np.loadtxt("../data/sim_snr.txt")
    t, res = np.loadtxt("../data/sim_ptr.txt", unpack=True)
    signal = np.loadtxt("../data/sim_signal.txt", unpack=True, usecols=1)

    Np = pulsars.shape[0]
    Nt = len(t) // Np
    signal = signal.reshape((Np, Nt))

    sd = pulsars[:, 4]
    resplot(t, res, sd, Np, Nt, signal, 0, 0, 'residual_sim.pdf')
    skymap(pulsars, sources, snr, None, 'skymap_sim.pdf')

    # fig, ax = plt.subplots()
    # ax.scatter(snr, 10**omega)
    # ax.axvline(30)
    # ax.axhline(np.pi / (14 / 365.25))
    # ax.axhline(np.pi / (14 * 130 / 365.25))
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.set_xlabel("SNR")
    # ax.set_ylabel("$\omega$")
    # plt.tight_layout()
    # fig.savefig("snr.pdf")
