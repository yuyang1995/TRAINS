import numpy as np
from param import parset
from plot_toolbox import par_names, identify, resplot, skymap, displot
from cluster import cluster
import sys

if __name__ == "__main__":
    sd = np.loadtxt("../{:s}".format(parset['PulsarCatalog']), usecols=4)
    t, res = np.loadtxt("../{:s}".format(parset['PTRFile']), unpack=True)
    res_con = np.loadtxt("../data/ptr_rec.txt")

    Np = len(sd)
    Nt = len(t) // Np
    Nps = len(res_con) // Np
    chi2 = np.zeros(Nps)
    for i in range(Nps):
        for j in range(Np):
            x = res_con[i * Np + j]
            y = res[j*Nt:j*Nt+Nt]
            chi2[i] += np.sum((x - y)**2 / sd[j]**2)
    chi2 /= (Np * Nt)
    ic_best = np.argmin(chi2)
    resplot(t, res, sd, Np, Nt, res_con, 0, ic_best)

    if parset['FlagMethod'] == '0':
        Ns = int(parset['SourceNumber'])
    else:
        Ns = 1
    sources = np.loadtxt("../data/sim_source.txt")
    snr = np.loadtxt("../data/sim_snr.txt")
    allsamples = np.loadtxt("../{:s}".format(parset['PosteriorSampleFile']),
                            usecols=range(8 * Ns))
    allsamples = allsamples.reshape((Nps, Ns, 8))

    if len(sys.argv) > 1:
        Nc = int(sys.argv[1])
        labels = cluster(allsamples, Nc)
        samples_list = [np.squeeze(allsamples[labels == i]) for i in range(Nc)]
    else:
        Nc = Ns
        samples_list = [allsamples[:, i, :] for i in range(Nc)]

    skymap(None, sources, snr, samples_list)

    mask = np.array(np.zeros(8 * Nc), dtype='bool')
    i = 0
    for c in parset['SourceParFix']:
        if i >= 8 * Ns:
            break
        if c == '1':
            mask[i] = True
        i = i + 1
    if parset['FlagEvolve'] == '0':
        for i in range(Nc):
            mask[i * 8 + 7] = True
    mask = mask.reshape((Nc, 8))

    for i in range(Nc):
        samples = samples_list[i]
        select = ~mask[i]
        index_truths = identify(samples, sources)
        truths = sources[index_truths]
        print(index_truths, snr[index_truths])

        displot(samples[:, select], par_names[select], truths[select],
                figname="corner_{:d}.pdf".format(i))
