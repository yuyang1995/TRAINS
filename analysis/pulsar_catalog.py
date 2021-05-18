import numpy as np
import sys
import scipy.io
from astropy.coordinates import SkyCoord

if __name__ == "__main__":
    Np = 100
    if len(sys.argv) > 1:
        Np = int(sys.argv[1])

    # read in the pulsar catalogue simulated for SKA
    # load input data
    skamsp = scipy.io.loadmat("survey_ska.mat")
    for key in skamsp:
        if not key.startswith('__'):
            skamsp[key] = skamsp[key].transpose()[0]

    # using Np randomly chosen pulsars
    # Np_tot = skamsp['D'].size
    # Ip = np.random.choice(Np_tot, Np)
    # using Np nearest pulsars
    Ip = np.argsort(skamsp['D'])[:Np]
    # sky location of the pulsars in the equatorial coordinate
    # we need to transfer from hr angle and degree to radian
    coord = SkyCoord(frame="galactic",
                     l=skamsp['l'][Ip], b=skamsp['b'][Ip], unit='deg')
    # right ascension and declination, in radian
    alphaP = coord.icrs.ra.rad
    deltaP = coord.icrs.dec.rad
    # (parallax) distance from SSB to pulsars, from kpc to ly
    distP = skamsp['D'][Ip] * 1e3 * 3.2615637771674333
    # 20% uncertainties for pulsar distances
    DeltadistP = distP * 0.2
    # 100 ns timing uncertainties
    sd = np.zeros(Np) + 1e-7
    # stack and save the pulsar information
    pulsar_info = np.stack((alphaP, deltaP, distP, DeltadistP, sd), axis=-1)
    np.savetxt("../data/pulsar_catalog.txt", pulsar_info, fmt=[
        '%.6f'] * 4 + ['%.6e'], header="{:d}".format(Np))
