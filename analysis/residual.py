import numpy as np
from scipy.integrate import quad


def FullResiduals(alphaS, sin_deltaS, cos_iota, psi, phi0, zeta, omega, tm,
                  alphaP, deltaP, t, distP, phiI, flag_sim=1, flag_evovle=0):
    '''
    Function to calculate the timing residuals for a full
    set of parameters (initial phase, arbitrary phase ...)
    '''
    # adjust the dimension to 3
    alphaS = alphaS[None, :, None]
    sin_deltaS = sin_deltaS[None, :, None]
    cos_iota = cos_iota[None, :, None]
    psi = psi[None, :, None]
    phi0 = phi0[None, :, None]
    zeta = zeta[None, :, None] * 0.5
    omega = omega[None, :, None]

    if flag_evovle:
        tm = tm[None, :, None]

    alphaP = alphaP[:, None, None]
    deltaP = deltaP[:, None, None]
    t = t[:, None, :]

    # antenna pattern functions (Np * Ns * 1)
    cos_deltaS = (1 - sin_deltaS**2)**0.5
    alphatilde = alphaS - alphaP
    Pp = -np.cos(deltaP)**2 * (1 - 2*np.cos(alphatilde)**2
                               + np.cos(alphatilde)**2*cos_deltaS**2) \
        + np.sin(deltaP)**2 * cos_deltaS**2 \
        - 0.5*np.sin(2*deltaP) * np.cos(alphatilde) * \
        2 * sin_deltaS * cos_deltaS
    Pc = -2*np.cos(deltaP)*np.sin(alphatilde) * \
        (np.cos(deltaP) * np.cos(alphatilde)*sin_deltaS
            - np.sin(deltaP) * cos_deltaS)
    cos_theta = np.cos(deltaP)*np.cos(alphaP)*cos_deltaS*np.cos(alphaS) + \
        np.cos(deltaP)*np.sin(alphaP)*cos_deltaS*np.sin(alphaS)\
        + np.sin(deltaP) * sin_deltaS
    Fp = Pp / (1 - cos_theta)
    Fc = Pc / (1 - cos_theta)

    # coefficient for each source (1 * Ns * 1)
    a1 = zeta * (1+cos_iota**2) * np.cos(2*psi)
    a2 = -zeta * 2*cos_iota * np.sin(2*psi)
    a3 = zeta * (1+cos_iota**2) * np.sin(2*psi)
    a4 = zeta * 2*cos_iota * np.cos(2*psi)

    # phase of pulsars (Np * Ns * 1)
    if flag_sim:
        distP = distP[:, None, None]
        tau = distP * (1 - cos_theta)
        if flag_evovle:
            phiI = phi0 + 0.8 * omega * tm * (1 - (1 + tau / tm)**0.625)
        else:
            phiI = phi0 - 0.5 * omega * tau
        phiI = np.fmod(phiI, np.pi)
    else:
        phiI = phiI[:, :, None]
        if flag_evovle:
            distP = distP[:, None, None]
            tau = distP * (1 - cos_theta)

    if flag_evovle:
        # Fourier coefficient
        sC = Fp * a2 + Fc * a4
        sS = Fp * a1 + Fc * a3

        # earth term
        beta = 1 - t / tm
        beta[beta < 0] = 0
        omegat = 2 * phi0 + 1.6 * omega * tm * (1 - beta**0.625)
        omegat = 2 * phi0 + omega * t
        res = (sC * np.cos(omegat) + sS * np.sin(omegat)) * beta**0.125

        # pulsar term
        omegap = omega * (1 + tau / tm)**-0.375
        tmp = tm + tau
        beta = 1 - t / tmp
        beta[beta < 0] = 0
        omegat = 2 * phiI + 1.6 * omegap * tmp * (1 - beta**0.625)
        res -= (sC * np.cos(omegat) + sS * np.sin(omegat)) * \
            beta**0.125 * (1 + tau / tm)**0.125
    else:
        DeltaC = np.cos(2*phi0)-np.cos(2*phiI)
        DeltaS = np.sin(2*phi0)-np.sin(2*phiI)

        sC = Fp * (DeltaC * a2 + DeltaS * a1) \
            + Fc * (DeltaC * a4 + DeltaS * a3)
        sS = Fp * (-DeltaS * a2 + DeltaC * a1) \
            + Fc * (-DeltaS * a4 + DeltaC * a3)

        omegat = omega * t
        res = sC * np.cos(omegat) + sS * np.sin(omegat)

    if flag_sim:
        res_mean = res.mean(axis=2, keepdims=True)
        res = res - res_mean
    else:
        res = np.einsum("ijk->ik", res, optimize=True)
        res_mean = res.mean(axis=1, keepdims=True)
        res = res - res_mean

    return res


def logLambda(phi, b):
    '''
    log likelihood ratio for a single pulsar
    '''
    value = b[0] + b[1] * np.cos(phi) + b[2] * np.sin(phi) + \
        b[3] * np.sin(2.0 * phi) + b[4] * np.cos(phi)**2
    return value


def Lambda(phi, b, b5):
    value = np.exp(logLambda(phi, b) - b5)
    return value


def LLR_Mx_Av(alphaS, sin_deltaS, cos_iota, psi, phi0, zeta, omega,
              alphaP, deltaP, sigmaP, t, res_data, method):
    Np = alphaP.size
    zeta = zeta * 0.5
    cos_deltaS = (1 - sin_deltaS**2)**0.5
    alphatilde = alphaS - alphaP
    Pp = -np.cos(deltaP)**2 * (1 - 2*np.cos(alphatilde)**2
                               + np.cos(alphatilde)**2*cos_deltaS**2) \
        + np.sin(deltaP)**2 * cos_deltaS**2 \
        - 0.5*np.sin(2*deltaP) * np.cos(alphatilde) * \
        2 * sin_deltaS * cos_deltaS
    Pc = -2*np.cos(deltaP)*np.sin(alphatilde) * \
        (np.cos(deltaP) * np.cos(alphatilde)*sin_deltaS
            - np.sin(deltaP) * cos_deltaS)
    cos_theta = np.cos(deltaP)*np.cos(alphaP)*cos_deltaS*np.cos(alphaS) + \
        np.cos(deltaP)*np.sin(alphaP)*cos_deltaS*np.sin(alphaS)\
        + np.sin(deltaP) * sin_deltaS
    Fp = Pp / (1 - cos_theta)
    Fc = Pc / (1 - cos_theta)

    A1 = (1 + cos_iota * cos_iota) * \
        (Fp * np.cos(2 * psi) + Fc * np.sin(2 * psi))
    A2 = 2 * cos_iota * (Fp * np.sin(2 * psi) - Fc * np.cos(2 * psi))
    A = 2 * zeta * np.sqrt(A1 * A1 + A2 * A2)
    phiA = np.arctan2(A1, A2)

    A = A[:, None]
    phiA = phiA[:, None]
    omegat = omega * t
    X = 0.5 * A * np.cos(phiA + omegat)
    Y = -0.5 * A * np.sin(phiA + omegat)
    Z = -0.5 * A * np.cos(2 * phi0 + phiA + omegat)

    sx = np.sum(res_data * X, axis=-1) / sigmaP**2
    sy = np.sum(res_data * Y, axis=-1) / sigmaP**2
    sz = np.sum(res_data * Z, axis=-1) / sigmaP**2
    xx = np.sum(X * X, axis=-1) / sigmaP**2
    xy = np.sum(X * Y, axis=-1) / sigmaP**2
    xz = np.sum(X * Z, axis=-1) / sigmaP**2
    yy = np.sum(Y * Y, axis=-1) / sigmaP**2
    yz = np.sum(Y * Z, axis=-1) / sigmaP**2
    zz = np.sum(Z * Z, axis=-1) / sigmaP**2

    c0 = -sx + xz
    c1 = sy - yz
    c2 = 0.5 * (xx - yy)
    c3 = -xy

    e4 = 4 * (c2 * c2 + c3 * c3)
    e3 = 4 * (c0 * c2 + c1 * c3)
    e2 = c0 * c0 + c1 * c1 - e4
    e1 = -2.0 * c1 * c3 - 4.0 * c0 * c2
    e0 = c3 * c3 - c0 * c0
    p = np.array([e0, e1, e2, e3, e4])

    b0 = sz - 0.5 * (yy + zz)
    b1 = sx - xz
    b2 = sy - yz
    b3 = -0.5 * xy
    b4 = 0.5 * (yy - xx)
    b = np.array([b0, b1, b2, b3, b4])

    roots = [np.roots(p[:, i]) for i in range(Np)]

    # calculate likehood and phiI at each root and boundary
    lh10 = np.zeros((10, Np)) - np.inf
    phiI10 = np.zeros((10, Np)) + np.nan
    # four roots
    for i in range(4):
        # y = cos(2 * phiI) must be real and its norm must <= 1
        rvalid = np.logical_and(np.imag(roots[i]) == 0, np.abs(roots[i]) <= 1)
        if np.any(rvalid):
            phiI10[i*2][rvalid] = np.arccos(roots[i][rvalid])
            lh10[i*2][rvalid] = logLambda(phiI10[i*2][rvalid], b[:, rvalid])
            phiI10[i*2+1][rvalid] = -phiI10[i*2][rvalid]
            lh10[i*2][rvalid] = logLambda(phiI10[i*2+1][rvalid], b[:, rvalid])
    # boundary
    phiI10[8] = np.zeros(Np)
    lh10[8] = logLambda(phiI10[8], b)
    phiI10[9] = np.zeros(Np) + np.pi
    lh10[9] = logLambda(phiI10[9], b)

    i0 = lh10.argmax(axis=0)
    i1 = np.indices((Np,))
    phiI = phiI10[i0, i1] / 2.0

    # choose the maximum likehood
    if method == 1:
        lh = lh10[i0, i1]
    elif method == 2:
        lh = np.empty(Np)
        b5 = lh10.max(axis=0)
        for i in range(Np):
            lh[i] = quad(Lambda, 0, 2 * np.pi, args=(b[:, i], b5[i]))[0]
    return lh, phiI
