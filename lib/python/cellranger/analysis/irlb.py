"""Copyright 2013 Michael Kane and Bryan Lewis.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from __future__ import annotations

# Based on the Python impl of IRLB by Bryan Lewis (https://github.com/bwlewis/irlbpy),
# but with support for sparse centering / scaling
import warnings

import numpy as np
import scipy.sparse as sp


# Matrix-vector product wrapper
# A is a numpy 2d array or matrix, or a scipy matrix or sparse matrix.
# x is a numpy vector only.
# Compute A.dot(x) if t is False,  A.transpose().dot(x)  otherwise.
def mult(A, x, t=False):
    assert x.ndim == 1
    if sp.issparse(A):
        if t:
            return sp.csr_matrix(x).dot(A).transpose().toarray()[:, 0]
        return A.dot(sp.csr_matrix(x).transpose()).toarray()[:, 0]
    if t:
        return np.asarray(A.transpose().dot(x)).ravel()
    return np.asarray(A.dot(x)).ravel()


def orthog(Y, X):
    """Orthogonalize a vector or matrix Y against the columns of the matrix X.

    This function requires that the column dimension of Y is less than X and
    that Y and X have the same number of rows.
    """
    dotY = Y.dot(X)
    return Y - mult(X, dotY)


# Simple utility function used to check linear dependencies during computation:
def invcheck(x):
    eps2 = 2 * np.finfo(float).eps
    if x > eps2:
        return 1.0 / x
    else:
        warnings.warn("Ill-conditioning encountered, result accuracy may be poor")
        return 0.0


def irlb(
    A: np.ndarray,
    n: int,
    tol: float = 0.0001,
    maxit: int = 50,
    center: np.ndarray | None = None,
    scale: np.ndarray | None = None,
    random_state=0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Estimate a few of the largest singular values and corresponding singular.

    vectors of matrix using the implicitly restarted Lanczos bidiagonalization
    method of Baglama and Reichel, see:

    Augmented Implicitly Restarted Lanczos Bidiagonalization Methods,
    J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005

    Keyword Arguments:
    tol   -- An estimation tolerance. Smaller means more accurate estimates.
    maxit -- Maximum number of Lanczos iterations allowed.

    Given an input matrix A of dimension j * k, and an input desired number
    of singular values n, the function returns a tuple X with five entries:

    X[0] A j * nu matrix of estimated left singular vectors.
    X[1] A vector of length nu of estimated singular values.
    X[2] A k * nu matrix of estimated right singular vectors.
    X[3] The number of Lanczos iterations run.
    X[4] The number of matrix-vector products run.

    The algorithm estimates the truncated singular value decomposition:
    A.dot(X[2]) = X[0]*X[1].
    """
    rs = np.random.RandomState(random_state)

    # Numpy routines do undesirable things if these come in as N x 1 matrices instead of size N
    # arrays
    if center is not None and not isinstance(center, np.ndarray):
        raise TypeError("center must be a numpy.ndarray")
    if scale is not None and not isinstance(scale, np.ndarray):
        raise TypeError("scale must be a numpy.ndarray")

    nu = n
    m = A.shape[0]
    n = A.shape[1]
    if min(m, n) < 2:
        raise ValueError("The input matrix must be at least 2x2.")
    # TODO: More efficient to have a path that performs a standard SVD
    # if over half the eigenvectors are requested
    m_b = min((nu + 20, 3 * nu, min(A.shape)))  # Working dimension size
    # m_b = nu + 7 # uncomment this line to check for similar results with R package
    mprod = 0
    it = 0
    j = 0
    k = nu
    smax = 1

    V: np.ndarray = np.zeros((n, m_b))  # Approximate right vectors
    W: np.ndarray = np.zeros((m, m_b))  # Approximate left vectors
    F: np.ndarray = np.zeros((n, 1))  # Residual vector
    B: np.ndarray = np.zeros((m_b, m_b))  # Bidiagonal approximation

    V[:, 0] = rs.randn(n)  # Initial vector
    V[:, 0] = V[:, 0] / np.linalg.norm(V)
    S: tuple[np.ndarray, np.ndarray, np.ndarray]

    while it < maxit:
        if it > 0:
            j = k

        VJ = V[:, j]

        # apply scaling
        if scale is not None:
            VJ = VJ / scale

        W[:, j] = mult(A, VJ)
        mprod += 1

        # apply centering
        # R code: W[, j_w] <- W[, j_w] - ds * drop(cross(dv, VJ)) * du
        if center is not None:
            W[:, j] = W[:, j] - np.dot(center, VJ)

        if it > 0:
            # NB W[:,0:j] selects columns 0,1,...,j-1
            W[:, j] = orthog(W[:, j], W[:, 0:j])
        s = np.linalg.norm(W[:, j])
        sinv = invcheck(s)
        W[:, j] = sinv * W[:, j]

        fn: float
        # Lanczos process
        while j < m_b:
            F = mult(A, W[:, j], t=True)
            mprod += 1

            # apply scaling
            if scale is not None:
                F = F / scale
            # apply centering, note for cases where center is the column
            # mean, this correction is often equivalent to a no-op as
            # np.sum(W[:, j]) is often close to zero
            if center is not None:
                sub = np.sum(W[:, j]) * center
                # This was in the R code, but not in python, allows test_pbmc4k_tiny to pass.
                if scale is not None:
                    sub = sub / scale
                F = F - sub
            F = F - s * V[:, j]
            F = orthog(F, V[:, 0 : (j + 1)])
            fn = np.linalg.norm(F)
            fninv = invcheck(fn)
            F = fninv * F
            if j < m_b - 1:
                V[:, j + 1] = F
                B[j, j] = s
                B[j, j + 1] = fn
                VJp1 = V[:, j + 1]

                # apply scaling
                if scale is not None:
                    VJp1 = VJp1 / scale

                W[:, j + 1] = mult(A, VJp1)
                mprod += 1

                # apply centering
                # R code: W[, jp1_w] <- W[, jp1_w] - ds * drop(cross(dv, VJP1)) * du
                if center is not None:
                    W[:, j + 1] = W[:, j + 1] - np.dot(center, VJp1)

                # One step of classical Gram-Schmidt...
                W[:, j + 1] = W[:, j + 1] - fn * W[:, j]
                # ...with full reorthogonalization
                W[:, j + 1] = orthog(W[:, j + 1], W[:, 0 : (j + 1)])
                s = np.linalg.norm(W[:, j + 1])
                sinv = invcheck(s)
                W[:, j + 1] = sinv * W[:, j + 1]
            else:
                B[j, j] = s
            j += 1
        # End of Lanczos process
        S = np.linalg.svd(B)
        R: np.ndarray = fn * S[0][m_b - 1, :]  # Residuals
        if it < 1:
            smax = S[1][0]  # Largest Ritz value
        else:
            smax = max((S[1][0], smax))

        conv = sum(np.abs(R[0:nu]) < tol * smax)
        if conv < nu:  # Not coverged yet
            k = max(conv + nu, k)
            k = min(k, m_b - 3)
        else:
            break
        # Update the Ritz vectors
        V[:, 0:k] = V[:, 0:m_b].dot(S[2].transpose()[:, 0:k])
        V[:, k] = F
        B = np.zeros((m_b, m_b))
        # Improve this! There must be better way to assign diagonal...
        for l in range(k):
            B[l, l] = S[1][l]
        B[0:k, k] = R[0:k]
        # Update the left approximate singular vectors
        W[:, 0:k] = W[:, 0:m_b].dot(S[0][:, 0:k])
        it += 1

    U: np.ndarray = W[:, 0:m_b].dot(S[0][:, 0:nu])
    V = V[:, 0:m_b].dot(S[2].transpose()[:, 0:nu])
    return (U, S[1][0:nu], V, it, mprod)
