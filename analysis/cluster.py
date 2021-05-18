import numpy as np
from sklearn.cluster import KMeans
# from sklearn.cluster import DBSCAN


def cluster(samples, Nc):
    Nps = samples.shape[0]
    Ns = samples.shape[1]
    samples = samples.reshape((Nps * Ns, 8))
    points = np.empty((Nps * Ns, 4))
    points[:, 0] = (1 - samples[:, 1]**2)**0.5 * np.cos(samples[:, 0])
    points[:, 1] = (1 - samples[:, 1]**2)**0.5 * np.sin(samples[:, 0])
    points[:, 2] = samples[:, 1]
    points[:, 3] = samples[:, 6] * 10  # omega

    clustering = KMeans(n_clusters=Nc, random_state=0).fit(points)
    # clustering = DBSCAN(eps=0.1, min_samples=Nc).fit(samples)
    labels = clustering.labels_

    return labels.reshape((Nps, Ns))
