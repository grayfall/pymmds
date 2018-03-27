"""

Classic metric MDS transforms a symmetric matrix of metric distances between N
objects into a set of projections in an orthonormal space. The method is akin to
PCA, though unlike PCA it doesn't provide a straightforward way to project new
objects onto the space. This package implements an active metric MDS method
described in "Principal Component and Correspondence Analysis in R" by
Dr. Herve Abdi (2017) and is thus similar to package bios2mds written in the R
language. The method allows new objects to be projected onto an existing space
defined by a set of active objects.

"""

from typing import cast

import numpy as np
import pandas as pd


__all__ = ["Space", 'read_dm']


class Space:
    def __init__(self, dm: pd.DataFrame):
        """
        Initialise the Space and calculate active samples' positions.
        The DM must be strictly finite, real and symmetric. It's index must be
        equal to its keys.
        :param dm: a metric distance matrix
        """
        if not len(dm):
            raise ValueError("a distance matrix can't be empty")
        if not (dm.keys() == dm.index).all():
            raise ValueError('row and column names do not match')
        if not len(set(dm.keys())):
            raise ValueError('duplicate object names')
        if len(dm.select_dtypes(include=[np.number]).columns) != len(dm.keys()):
            raise ValueError('a distance matrix must be strictly numeric')
        distances = cast(np.ndarray, dm.as_matrix().copy())
        if not np.isfinite(distances).all():
            raise ValueError("all values in the distance matrix must be finite")
        if distances.diagonal().any():
            raise ValueError("a distance matrix can't have nonzero diagonal "
                             "elements")
        if not (distances == distances.T).all():
            raise ValueError("a distance matrix must be absolutely symmetric")
        # calculate coordinates for active samples
        n = len(dm.keys())
        d = distances ** 2
        masses = np.identity(n) - np.full([n]*2, 1/n)
        s = -0.5 * (masses @ d @ masses.T)
        # decompose
        eigen = np.linalg.eigh(s)
        # drop negative eigenvalues and sort by eigenvalues in descending order
        ordering = np.where(eigen[0] > 0)[0][::-1]
        values = eigen[0][ordering]
        vectors = eigen[1][:, ordering]
        # compute and format coordinates
        coord = np.diag(np.repeat((1/n)**(-0.5), n)) @ vectors @ np.diag(values**0.5)
        # keep items required to project supplementary values
        self._keys = list(dm.keys())
        self._active = pd.DataFrame(coord, index=dm.keys())
        self._d_act = d
        self._masses_act = masses
        self._values = values
        self._vectors = vectors
        self._explained = (values * 100 / values.sum()).round(3)

    @property
    def explained(self) -> np.ndarray:
        """
        The fraction of variance explained by each dimension.
        :return: an array of floating point numbers
        """
        return self._explained.copy()

    @property
    def keys(self) -> list:
        """
        Active samples
        :return: a list of names
        """
        return self._keys[:]

    @property
    def ndim(self) -> int:
        """
        The number of dimensions. There can be at most `len(keys)` dimensions,
        though it is not uncommon to get fewer dimensions. For example, due to
        numerical instability or poor metric estimations some dimensions might
        end up with near-zero negative eigenvalues and are consequently dropped.
        :return: an integer
        """
        return self._active.shape[1]

    @property
    def active(self) -> pd.DataFrame:
        """
        Active samples's coordinates.
        :return: a pandas data frame with coordinates of active samples,
        one sample per row; each column encodes a dimension in the Space;
        columns are sorted with respect to the fraction of variance explained
        by the corresponding dimensions in descending order
        """
        return self._active.copy()

    def project(self, dm: pd.DataFrame) -> pd.DataFrame:
        """
        Project supplementary samples onto the Space. The DM must be strictly
        finite and real. There should be a column for each active sample.
        :param dm: a metric distance matrix
        :return: a pandas data frame with coordinates of supplementary samples,
        one sample per row; each column encodes a dimension in the Space;
        columns are sorted with respect to the fraction of variance explained
        by the corresponding dimensions in descending order.
        """
        # make sure `dm` is fine
        if not len(dm):
            raise ValueError("a distance matrix can't be empty")
        if not set(dm.keys()) == set(self.keys):
            raise ValueError('a dm must contain distances from supplementary '
                             'samples to active samples and its columns must '
                             'be named after the active samples')
        if len(dm.select_dtypes(include=[np.number]).columns) != len(dm.keys()):
            raise ValueError('a distance matrix must be strictly numeric')
        # make sure all columns are in correct order
        distances = cast(np.ndarray, dm[self.keys].as_matrix().copy())
        if not np.isfinite(distances).all():
            raise ValueError("all values in the distance matrix must be finite")
        n_act = len(self.keys)
        n_sup = distances.shape[0]
        d_sup = distances ** 2
        masses_sup = np.full((n_act, n_sup), (1 / n_act))
        s_sup = -0.5 * self._masses_act @ (d_sup.T - (self._d_act @ masses_sup))
        f_sup = s_sup.T @ self.active.as_matrix() @ np.diag(self._values**-1)
        return pd.DataFrame(f_sup, index=list(dm.index))


def read_dm(path: str, sep: str='\t') -> pd.DataFrame:
    """
    Read a distance matrix. The first row and column are used as column and row
    indices respectively. The function can parse both active and supplementary
    DMs.
    :param path: a path to a (possible gzipped) table; compressing is inferred
    from file extension.
    :param sep: table separator symbol
    :return: a pandas data frame (dm); dm.keys() and dm.index contain sample
    names
    """
    # a pandas bug renders `dtype`/`converters` useless for `index_col`, making
    # it necessary to create an intermediate DF
    raw = pd.read_csv(path, sep=sep, dtype={0: str})
    return raw.set_index(raw.keys()[0])


if __name__ == '__main__':
    raise RuntimeError
