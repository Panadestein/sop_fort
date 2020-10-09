"""
The present module decomposes a reference tensor in Tucker form
using the Tensorly library
"""
import time
from ctypes import CDLL, c_double
import numpy as np
from numpy.polynomial import chebyshev as cheby
import tensorly as tl
from rho_sop import soputils

# Regularize coulomb (C Function)

CFUNC = CDLL("./coul.so")
reg_coulomb = CFUNC.reg_coul
reg_coulomb.restype = c_double
reg_coulomb.argtypes = [c_double, c_double, c_double, c_double,
                        c_double, c_double]

# Generate ab initio reference grid

X1MIN, X1MAX = -4.5, 4.5
Y1MIN, Y1MAX = -4.5, 4.5
Z1MIN, Z1MAX = -4.5, 4.5
X2MIN, X2MAX = -4.5, 4.5
Y2MIN, Y2MAX = -4.5, 4.5
Z2MIN, Z2MAX = -4.5, 4.5

X1 = np.linspace(X1MIN, X1MAX, num=2)
Y1 = np.linspace(Y1MIN, Y1MAX, num=2)
Z1 = np.linspace(Z1MIN, Z1MAX, num=2)
X2 = np.linspace(X2MIN, X2MAX, num=2)
Y2 = np.linspace(Y2MIN, Y2MAX, num=2)
Z2 = np.linspace(Z2MIN, Z2MAX, num=2)

X1_M, Y1_M, Z1_M, X2_M, Y2_M, Z2_M = np.meshgrid(
    X1, Y1, Z1, X2, Y2, Z2, indexing='ij')

E_AB = np.vectorize(reg_coulomb)(X1_M, Y1_M, Z1_M, X2_M, Y2_M, Z2_M)

GRID = np.meshgrid(X1, Y1, Z1, X2, Y2, Z2, indexing='ij')
INTCOORD = np.vstack(list(map(np.ravel, GRID))).T

# SOP-FBR function

CHEBDIM = np.array([10, 10, 10, 10, 10, 10])
CTEN_DIM = np.array([5, 5, 5, 5, 5, 5])
NCHEB = np.sum(CTEN_DIM * CHEBDIM)
CARRAY = np.loadtxt('params_sop')
CHEB = CARRAY[:NCHEB]
CORA = CARRAY[NCHEB:].reshape(CTEN_DIM)


def vchpot(*q_array):
    """Computes the SOP potential for the reference geometries using the
    tensor n-mode approach"""

    v_matrices = []
    idx_cheb = 0
    for kdof, m_kp in enumerate(CTEN_DIM):
        v_kp = np.zeros(m_kp)
        for j_kp in np.arange(m_kp):
            v_kp[j_kp] = cheby.chebval(q_array[kdof],
                                       CHEB[idx_cheb:idx_cheb + CHEBDIM[kdof]])
            idx_cheb += CHEBDIM[kdof]
        v_matrices.append(v_kp)

    prod = tl.tucker_tensor.tucker_to_tensor((CORA, v_matrices))[0]

    return prod


# Evaluate potential and compute RMSE

T0 = time.time()
E_SOP = np.vectorize(vchpot)(X1_M, Y1_M, Z1_M, X2_M, Y2_M, Z2_M)
T1 = time.time()
RMS = np.sqrt(((E_AB - E_SOP) ** 2).mean())
print(f"RMSE {RMS}")
print(f"Time {T1 - T0}")

# Try F2PY wrapper

T0F = time.time()
RMSF = soputils.rho(CARRAY, INTCOORD, E_AB.flatten(), CTEN_DIM, CHEBDIM)
RMSF_CHEB = soputils.rho_cheb(CHEB, CORA.flatten(), INTCOORD,
                              E_AB.flatten(), CTEN_DIM, CHEBDIM)
RMSF_CORE = soputils.rho_core(CORA.flatten(), CHEB, INTCOORD,
                              E_AB.flatten(), CTEN_DIM, CHEBDIM)
T1F = time.time()
print(f"RMSE {RMSF}")
print(f"RMSE_CHEB {RMSF_CHEB}")
print(f"RMSE_CORE {RMSF_CORE}")
print(f"Time {T1F - T0F}")
