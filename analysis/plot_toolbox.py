import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import corner
import astropy.units as units
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord


plt.rc('font', family='serif')
mpl.rcParams['text.usetex'] = 'True'
mpl.rcParams['font.size'] = 15
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

par_names = [r'$\alpha$', r'$\sin \delta$', r'$\cos \iota$',
             r'$\psi$', r'$\phi$', r'$\log \zeta$',
             r'$\log \omega$', r'$\log t_{\rm merge}$']
par_names = np.array(par_names)


def identify(samples, sources, indexes=[0, 1, 5, 6]):
    # alpha, sin(delta), log(zeta) and omega are used
    # to identify the source defautly

    chi2 = np.zeros(sources.shape[0])
    for i in indexes:
        mean = np.mean(samples[:, i])
        std = np.std(samples[:, i])
        diff = np.abs(sources[:, i] - mean)
        # note the topology of alpha
        if i == 0:
            diff = np.minimum(diff, 2 * np.pi - diff)
        if std == 0:
            std = 1e-10
        chi2 += (diff / std)**2

    return chi2.argmin()


def resplot(t, res, sd, Np, Nt, res_con=None, Nps=0,
            ic_best=None, figname='residuals.pdf'):
    with PdfPages(figname) as pdf:
        for i in range(Np):
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.errorbar(t[i*Nt:i*Nt+Nt], res[i*Nt:i*Nt+Nt],
                        yerr=sd[i], fmt='o')
            if res_con is not None:
                for j in range(Nps):
                    ax.plot(t[i*Nt:i*Nt+Nt], res_con[j*Np+i],
                            lw=1, c='k', alpha=0.2)
                if ic_best is not None:
                    ax.plot(t[i*Nt:i*Nt+Nt], res_con[ic_best*Np+i],
                            lw=2, c='r')

            ax.set_xlabel(r'time (yr)')
            ax.set_ylabel(r'residual(s)')
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def skymap(pulsars, sources, snr, samples_list, figname='skymap.pdf'):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection="mollweide")

    if pulsars is not None:
        galactic_plane = SkyCoord(frame="galactic", l=np.linspace(
            -63, 296, 100), b=np.zeros(100), unit='deg')
        alpha = galactic_plane.icrs.ra.wrap_at(180 * units.degree).radian
        delta = galactic_plane.icrs.dec.rad
        ax.plot(alpha, delta, c='k', alpha=0.5, ls='--', lw=0.8)

        alpha = coord.Angle(
            pulsars[:, 0] * units.rad).wrap_at(180 * units.degree).radian
        ax.scatter(alpha, pulsars[:, 1], marker='^', c='k', s=5)

    if samples_list is not None:
        Nc = len(samples_list)
        log_omega = np.empty(Nc)
        for i in range(Nc):
            samples = samples_list[i]
            log_omega[i] = np.mean(samples[:, 6])
        index = np.argsort(log_omega)

        for j in range(Nc):
            i = index[j]
            if Nc > 1:
                color = mpl.cm.gist_rainbow(j / (Nc - 1))
            else:
                color = mpl.cm.gist_rainbow(j)
            samples = samples_list[i]
            alpha = coord.Angle(
                samples[:, 0] * units.rad).\
                wrap_at(180 * units.degree).radian
            delta = np.arcsin(samples[:, 1])
            ax.scatter(alpha, delta, s=1,
                       label='{:d}'.format(i), color=color)

    if sources is not None:
        alpha = coord.Angle(
            sources[:, 0] * units.rad).wrap_at(180 * units.degree).radian
        delta = np.arcsin(sources[:, 1])
        ax.scatter(alpha, delta, c=np.log(snr), s=20, cmap='Blues')

    ax.set_xticklabels(['14h', '16h', '18h', '20h', '22h',
                        '0h', '2h', '4h', '6h', '8h', '10h'])

    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')
    ax.grid(True)
    plt.tight_layout()
    fig.savefig(figname)


def displot(samples, labels, truths=None, estimates=None,
            figname='corner.pdf'):
    levels = 1.0 - np.exp(-0.5 * np.linspace(1.0, 2.0, 3)**2)
    sample_max = np.max(samples, axis=0)
    sample_min = np.min(samples, axis=0)
    ranges = np.stack((sample_min, sample_max), axis=-1)

    fig = corner.corner(samples, quantiles=[0.16, 0.5, 0.84], show_titles=True,
                        title_kwargs={"fontsize": 15}, labels=labels,
                        label_kwargs={"fontsize": 18}, smooth=True,
                        range=ranges, levels=levels)

    if truths is not None:
        corner.overplot_lines(fig, truths, color="g")
        corner.overplot_points(fig, truths[None], marker="s", color="g")
    if estimates is not None:
        if estimates == 'mean':
            estimates = np.mean(samples, axis=0)
        elif isinstance(estimates, int):
            estimates = samples[estimates]

        corner.overplot_lines(fig, estimates, color="r")
        corner.overplot_points(fig, estimates[None], marker="s", color="r")
    fig.savefig(figname)
