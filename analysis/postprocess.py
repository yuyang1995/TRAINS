import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import getopt


def logsumexp(values):     # log SUM( exp(values) )
    biggest = np.max(values)
    x = values - biggest
    result = np.log(np.sum(np.exp(x))) + biggest
    return result


def logdiffexp(x1, x2):    # log( exp(x1) - exp(x2))
    biggest = x1
    xx1 = x1 - biggest
    xx2 = x2 - biggest
    result = np.log(np.exp(xx1) - np.exp(xx2)) + biggest
    return result


def postprocess():
    cut = 0
    numResampleLogX = 1
    compression_bias_min = 1.
    compression_scatter = 0.
    zoom_in = True
    moreSamples = 1

    str_mode = '_pt'
    temperature = 0
    opts, args = getopt.getopt(sys.argv[1:], "t:l:")
    for op, value in opts:
        if op == "-t":
            temperature = float(value)
        elif op == "-l":
            str_mode = value

    pdf = PdfPages('dnest' + str_mode + '.pdf')

    levels_orig = np.genfromtxt(
        "../data/levels"+str_mode+".txt", comments='#',
        dtype=float, skip_footer=1)
    sample_info = np.genfromtxt(
        "../data/sample_info"+str_mode+".txt", comments='#',
        dtype=float, skip_footer=1)
    sample = np.atleast_2d(np.genfromtxt(
        "../data/sample"+str_mode+".txt", comments='#',
        skip_footer=1, dtype=float))

    # reset level assignment
    # idx = (sample_info[:, 0] > levels_orig.shape[0] - 1)
    # sample_info[idx, 0] = levels_orig.shape[0] - 1

    # sample_info = sample_info[10000:,:]
    # sample = sample[10000:,:]

    print(levels_orig.shape, sample.shape, sample_info.shape)
    with open("../data/sampler_state"+str_mode+".txt", "w") as f:
        new_line = "{0:d}\t{1:d}\n".format(
            levels_orig.shape[0] + 1, sample_info.shape[0] + 1)
        f.write(new_line)

    idx = np.where(sample_info[:, 0] < levels_orig.shape[0])
    sample = sample[idx[0], :]
    sample_info = sample_info[idx[0], :]
    sample = sample[int(cut*sample.shape[0]):, :]
    sample_info = sample_info[int(cut*sample_info.shape[0]):, :]

    if sample.shape[0] != sample_info.shape[0]:
        print('# Size mismatch. Truncating...')
    lowest = np.min([sample.shape[0], sample_info.shape[0]])
    sample = sample[0:lowest, :]
    sample_info = sample_info[0:lowest, :]

    plt.figure(1)
    plt.plot(sample_info[:, 0], "k")
    plt.xlabel("Iteration")
    plt.ylabel("Level")

    pdf.savefig()

    plt.figure(2)
    plt.subplot(2, 1, 1)
    plt.plot(np.diff(levels_orig[:, 0]), "k")
    plt.ylabel("Compression")
    plt.xlabel("Level")
    xlim = plt.gca().get_xlim()
    plt.axhline(-1., color='g')
    plt.axhline(-np.log(10.), color='g', linestyle="--")
    plt.ylim(top=0.05)

    plt.subplot(2, 1, 2)
    good = np.nonzero(levels_orig[:, 4] > 0)[0]
    plt.plot(levels_orig[good, 3]/levels_orig[good, 4], marker='o')
    plt.xlim(xlim)
    plt.ylim([0., 1.])
    plt.xlabel("Level")
    plt.ylabel("MH Acceptance")
    pdf.savefig()

    logl_levels = [(levels_orig[i, 1], levels_orig[i, 2])
                   for i in range(0, levels_orig.shape[0])]
    logl_samples = [(sample_info[i, 1], sample_info[i, 2], i)
                    for i in range(0, sample.shape[0])]
    logx_samples = np.zeros((sample_info.shape[0], numResampleLogX))
    logp_samples = np.zeros((sample_info.shape[0], numResampleLogX))
    logP_samples = np.zeros((sample_info.shape[0], numResampleLogX))
    P_samples = np.zeros((sample_info.shape[0], numResampleLogX))
    logz_estimates = np.zeros((numResampleLogX, 1))
    H_estimates = np.zeros((numResampleLogX, 1))

    sandwich = sample_info[:, 0].copy().astype('int')
    for i in range(0, sample.shape[0]):
        while sandwich[i] < levels_orig.shape[0] - 1 and \
                logl_samples[i] > logl_levels[sandwich[i] + 1]:
            sandwich[i] += 1

    for z in range(0, numResampleLogX):
        levels = levels_orig.copy()
        compressions = -np.diff(levels[:, 0])
        compressions *= compression_bias_min + \
            (1.0 - compression_bias_min) * np.random.rand()
        compressions *= np.exp(compression_scatter *
                               np.random.randn(compressions.size))
        levels[1:, 0] = -compressions
        levels[:, 0] = np.cumsum(levels[:, 0])

        for i in range(0, levels.shape[0]):
            which = np.nonzero(sandwich == i)[0]
            logl_samples_thisLevel = []
            for j in range(0, len(which)):
                logl_samples_thisLevel.append(
                    copy.deepcopy(logl_samples[which[j]]))
            logl_samples_thisLevel = sorted(logl_samples_thisLevel)

            N = len(logl_samples_thisLevel)

            logx_max = levels[i, 0]
            if i == levels.shape[0] - 1:
                logx_min = -1E300   # the minmial X -> 0
            else:
                logx_min = levels[i+1, 0]

            Umin = np.exp(logx_min - logx_max)

            if N == 0 or numResampleLogX > 1:
                U = Umin + (1.0 - Umin) * np.random.rand(len(which))
            else:
                U = Umin + (1.0 - Umin) * np.linspace(
                    1.0 / (N + 1), 1.0 - 1.0/(N+1), N)

            logx_samples_thisLevel = np.sort(logx_max + np.log(U))[::-1]

            for j in range(0, which.size):
                logx_samples[logl_samples_thisLevel[j]
                             [2]][z] = logx_samples_thisLevel[j]

                if j != which.size - 1:
                    left = logx_samples_thisLevel[j+1]
                elif i == levels.shape[0] - 1:
                    left = -1E300
                else:
                    left = levels[i+1][0]

                if j != 0:
                    right = logx_samples_thisLevel[j-1]
                else:
                    right = levels[i][0]

                logp_samples[logl_samples_thisLevel[j][2]][z] = np.log(
                    0.5) + logdiffexp(right, left)

        logp_samples[:, z] = logp_samples[:, z] - \
            logsumexp(logp_samples[:, z])  # make sure that Sum(p) = 1.0

        if temperature == 0:
            temperature = -(levels_orig[-1][1] - levels_orig[-2]
                            [1]) / (levels_orig[-1][0] - levels_orig[-2][0])
            temperature = max(1.0, temperature * 1.2)
            print('The proper temperature is {:.2f}'.format(temperature))
        else:
            print('Set the temperature to {:.2f}'.format(temperature))

        logl = sample_info[:, 1]/temperature
        logP_samples[:, z] = logp_samples[:, z] + logl
        logz_estimates[z] = logsumexp(logP_samples[:, z])
        logP_samples[:, z] -= logz_estimates[z]
        P_samples[:, z] = np.exp(logP_samples[:, z])
        H_estimates[z] = -logz_estimates[z] + np.sum(P_samples[:, z]*logl)

        plt.figure(3)
        plt.subplot(2, 1, 1)

        plt.plot(logx_samples[:, z], sample_info[:, 1], 'b.', label='Samples')

        plt.plot(levels[1:, 0], levels[1:, 1], 'r.', label='Levels')
        plt.legend(numpoints=1, loc='lower left')
        plt.ylabel('log(L)')
        plt.title(str(z+1) + "/" + str(numResampleLogX) +
                  ", log(Z) = " + str(logz_estimates[z][0]))
        # Use all plotted logl values to set ylim
        combined_logl = np.hstack([sample_info[:, 1], levels[1:, 1]])
        combined_logl = np.sort(combined_logl)
        lower = combined_logl[int(0.1*combined_logl.size)]
        upper = combined_logl[-1]
        diff = upper - lower
        lower -= 0.05*diff
        upper += 0.05*diff
        if zoom_in:
            plt.ylim([lower, upper])

        xlim = plt.gca().get_xlim()

        plt.subplot(2, 1, 2)

        plt.plot(logx_samples[:, z], P_samples[:, z], 'b.')
        plt.ylabel('Posterior Weights')
        plt.xlabel('log(X)')
        plt.xlim(xlim)
        ylim = plt.gca().get_ylim()
        plt.text(xlim[0]+0.1*(xlim[1]-xlim[0]), ylim[1] -
                 0.1*(ylim[1]-ylim[0]), 'T='+str(temperature))
        pdf.savefig()

    pdf.close()
    P_samples = np.mean(P_samples, 1)
    P_samples = P_samples/np.sum(P_samples)
    logz_estimate = np.mean(logz_estimates)
    logz_error = np.std(logz_estimates)
    H_estimate = np.mean(H_estimates)
    H_error = np.std(H_estimates)
    ESS = np.exp(-np.sum(P_samples*np.log(P_samples+1E-300)))

    print("log(Z) = " + str(logz_estimate) + " +- " + str(logz_error))
    print("Information = " + str(H_estimate) +
          " +- " + str(H_error) + " nats.")
    print("Effective sample size = " + str(ESS))

    # plt.show()

    N = int(moreSamples*ESS)
    posterior_sample = np.zeros((N, sample.shape[1]))
    w = P_samples
    w = w/np.max(w)

    np.savetxt('../data/weigths'+str_mode+'.txt', w)

    for i in range(0, N):
        while True:
            which = np.random.randint(sample.shape[0])
            if np.random.rand() < w[which]:
                break
        posterior_sample[i, :] = sample[which, :]

    np.savetxt('../data/posterior_sample'+str_mode+'.txt', posterior_sample)
    return temperature


if __name__ == "__main__":
    postprocess()
