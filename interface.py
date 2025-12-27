import numpy as np
import torch
from scipy.sparse.linalg import LinearOperator, lsmr
    
def Create_ComMat(T, n):
    C = torch.zeros(T*n, T*n, dtype=torch.double)
    for i in range(T):
        for j in range(n):
            e_i = torch.zeros(T, 1)
            e_i[i] = 1
            e_j = torch.zeros(n, 1)
            e_j[j] = 1
            C += torch.kron(e_i @ e_j.T, e_j @ e_i.T)
    return C

def Jv(p, v, F):
    _, jvp = torch.autograd.functional.jvp(F, p, v)
    return jvp

def JTv(p, v, F):
    _, vjp = torch.autograd.functional.vjp(F, p, v)
    return vjp

def make_linear_operator_for_J(p, F):
    n = p.numel()
    m = F(p).numel()

    def matvec(x_np):
        x = torch.tensor(x_np, dtype=p.dtype)
        Jx = JTv(p, x, F)         # Note: vJ for transpose
        return Jx.detach().numpy()

    def rmatvec(y_np):
        y = torch.tensor(y_np, dtype=p.dtype)
        JTy = Jv(p, y, F)         # Note: Jv for transpose
        return JTy.detach().numpy()

    return LinearOperator(
        shape=(n, m),
        matvec=matvec,
        rmatvec=rmatvec,
        dtype=np.float64
    )

def solve_rectangular_with_LSMR(p, b, F, damp=0.00):
    A = make_linear_operator_for_J(p, F)
    b = b.detach().numpy()
    result = lsmr(A, b, damp=damp)
    x = result[0]  # solution
    return x

def WTdLdD(L, S, Lam, D, Q, R, gamma, n, m, T, V, Pi, WT):

    def _to_torch(x, force_2d=False):
        arr = np.array(x, dtype=np.float64, copy=True)
        if arr.ndim == 0:
            arr = arr.reshape(1, 1) if force_2d else arr.reshape(1)
        return torch.from_numpy(arr).to(torch.double)
    
    torch.set_default_dtype(torch.float64)
    
    n_i = int(n); m_i = int(m); T_i = int(T)
    barE1 = torch.cat([torch.eye(n_i+m_i, dtype=torch.double), torch.zeros(2*n_i, n_i+m_i, dtype=torch.double)], dim=0)
    barE2 = torch.cat([torch.zeros(n_i+m_i, 2*n_i, dtype=torch.double), torch.eye(2*n_i, dtype=torch.double)], dim=0)
    E1_t = torch.kron(barE1, barE1)
    E2_t = torch.kron(barE2, barE2)
    tmp = torch.cat([torch.eye(m_i, dtype=torch.double), torch.zeros(n_i, m_i, dtype=torch.double)], dim=0)
    dF1dS_t = torch.kron(tmp, tmp)
    dF2dS_t = torch.zeros(4*n_i**2, m_i**2, dtype=torch.double)

    S_t   = _to_torch(S, force_2d=True)
    L_t   = _to_torch(L, force_2d=True)
    Lam_t = _to_torch(Lam, force_2d=True)
    D_t   = _to_torch(D, force_2d=True)
    Z_t   = _to_torch(D_t[0:n_i, :], force_2d=True)
    X_t   = _to_torch(D_t[n_i:2*n_i, :], force_2d=True)
    U_t   = _to_torch(D_t[2*n_i:, :], force_2d=True)
    Q_t   = _to_torch(Q, force_2d=True)
    R_t   = _to_torch(R, force_2d=True)
    V_t   = _to_torch(V, force_2d=True)
    Pi_t  = _to_torch(Pi, force_2d=True)
    gamma_f = float(gamma)

    C = Create_ComMat(T_i, n_i)

    # --- Residual function (returns flattened G tensor) ---
    def residual(L_var, S_var, Lam_var, D_var):
        C = Create_ComMat(T_i, n_i)
        Z_var = D_var[0:n_i, :]
        X_var = D_var[n_i:2*n_i, :]
        U_var = D_var[2*n_i:, :]

        # dF1/dL
        t1 = torch.cat([torch.zeros(m_i, n_i), torch.eye(n_i)], dim=0)
        t2 = torch.cat([V_t @ U_var, torch.zeros(n_i, T_i)], dim=0)
        dF1dL_loc = (torch.kron(t1, t2) +
                     torch.kron(t2, t1) @ C +
                     torch.kron(t1, torch.cat([torch.zeros(m_i, T_i), X_var], dim=0)))
        # dF2/dL
        t1b = torch.cat([torch.eye(n_i, dtype=torch.double), torch.zeros(n_i, n_i, dtype=torch.double)], dim=0)
        t2b = torch.cat([torch.zeros(n_i, n_i, dtype=torch.double), torch.eye(n_i, dtype=torch.double)], dim=0)
        t3b = torch.cat([Z_var, torch.zeros(n_i, T_i, dtype=torch.double)], dim=0)
        dF2dL_loc = (torch.kron(t1b, torch.cat([X_var, torch.zeros(n_i, T_i)], dim=0)) +
                     torch.kron(t2b, t3b) +
                     torch.kron(t3b, t2b) @ C +
                     torch.kron(t2b, torch.cat([torch.zeros(n_i, T_i), X_var], dim=0)))
        dFdL_loc = E1_t @ dF1dL_loc + E2_t @ dF2dL_loc
        dFdS_loc = E1_t @ dF1dS_t + E2_t @ dF2dS_t

        #tmp1 = (Q_t @ X_var).T.flatten()
        #tmp2 = (Pi_t @ L_var).flatten()
        #tmp3 = Lam_var.flatten()
        tmp1 = (X_var.T @ Q_t.T).T.reshape(-1)
        tmp2 = (Pi_t @ L_var).T.reshape(-1)
        tmp3 = Lam_var.T.reshape(-1)
        G1_loc = tmp1 + 2 * gamma_f * tmp2 - tmp3 @ dFdL_loc
        tmpS = torch.eye(m_i, dtype=torch.double).flatten()
        G2_loc = tmpS - tmp3 @ dFdS_loc

        VUL = V_t @ U_var @ L_var
        XL = X_var @ L_var
        ZL = Z_var @ L_var
        F1_upper = torch.cat([S_var, VUL], dim=1)
        F1_lower = torch.cat([VUL.T, XL], dim=1)
        F1 = torch.cat([F1_upper, F1_lower], dim=0)
        F2_upper = torch.cat([XL - torch.eye(n_i, dtype=torch.double), ZL], dim=1)
        F2_lower = torch.cat([ZL.T, XL], dim=1)
        F2 = torch.cat([F2_upper, F2_lower], dim=0)
        zeros_12 = torch.zeros(m_i + n_i, 2 * n_i, dtype=torch.double)
        zeros_21 = torch.zeros(2 * n_i, m_i + n_i, dtype=torch.double)
        F_top = torch.cat([F1, zeros_12], dim=1)
        F_bottom = torch.cat([zeros_21, F2], dim=1)
        F_mat = torch.cat([F_top, F_bottom], dim=0)

        G3_loc = F_mat @ Lam_var.T
        G4_loc = Lam_var - Lam_var.T
        G5_loc = Lam_var[:m_i + n_i, m_i + n_i:]
        G6_loc = Lam_var[m_i + n_i:, :m_i + n_i]
        return torch.cat([
           G1_loc.T.reshape(-1), G2_loc.T.reshape(-1), G3_loc.T.reshape(-1),
           G4_loc.T.reshape(-1), G5_loc.T.reshape(-1), G6_loc.T.reshape(-1)
        ])

    def res_vec_LSLam(LSLam_vec, D_var, LT_s, ST_s, LamT_s):
        L_vec = LSLam_vec[0:LT_s[0]*LT_s[1]]
        S_vec = LSLam_vec[LT_s[0]*LT_s[1]:LT_s[0]*LT_s[1]+ST_s[0]*ST_s[1]]
        Lam_vec = LSLam_vec[LT_s[0]*LT_s[1]+ST_s[0]*ST_s[1]:]
        LT = L_vec.reshape(LT_s)
        ST = S_vec.reshape(ST_s)
        LamT = Lam_vec.reshape(LamT_s)
        res = residual(LT.T, ST.T, LamT.T, D_var)
        return res

    def res_vec_D(L_var, S_var, Lam_var, D_vec, DT_s):
        DT = D_vec.reshape(DT_s)
        res = residual(L_var, S_var, Lam_var, DT.T)
        return res

    L_vec = L_t.T.reshape(-1)
    LT_s = L_t.T.shape
    S_vec = S_t.T.reshape(-1)
    ST_s = S_t.T.shape
    Lam_vec = Lam_t.T.reshape(-1)
    LamT_s = Lam_t.T.shape
    LSLam_vec = torch.cat([L_vec, S_vec, Lam_vec])
    D_vec = D_t.T.reshape(-1)
    DT_s = D_t.T.shape
    zeros1 = torch.zeros(S_t.numel())
    zeros2 = torch.zeros(Lam_t.numel())
    WT_t = torch.from_numpy(WT).to(torch.double)
    rhs = torch.cat([WT_t, zeros1, zeros2])

    G_v = lambda v: res_vec_LSLam(v, D_t, LT_s, ST_s, LamT_s)
    xi=solve_rectangular_with_LSMR(LSLam_vec, rhs, G_v, damp=0.0)
    xi_t = torch.from_numpy(xi).to(torch.double)
    G_D = lambda d: res_vec_D(L_t, S_t, Lam_t, d, DT_s)
    _, WTdLdD_adj = torch.autograd.functional.vjp(G_D, D_vec, -xi_t)

    return WTdLdD_adj