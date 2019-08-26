import numpy as np
from scipy.special import binom
from scipy.sparse import issparse
from MKLpy.utils.validation import check_X_T
import types

def rbf_kernel(X, T=None, sigma=1):
    """performs the HPK kernel between the samples matricex *X* and *T*.
    The HPK kernel is defines as:
    .. math:: k(x,z) = \langle x,z \rangle^d

    Parameters
    ----------
    X : (n,m) array_like,
        the train samples matrix.
    T : (l,m) array_like,
        the test samples matrix. If it is not defined, then the kernel is calculated
        between *X* and *X*.

    Returns
    -------
    K : (l,n) ndarray,
        the HPK kernel matrix.
    """

    X, T = check_X_T(X, T)
    return np.exp(-(X - T)**2/(2*sigma**2))