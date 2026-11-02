# Degree-6 sector table (p = 10^9+7; sector key = character w.r.t. the engine generators (3,3,3,0,0,0), (5,2,0,1,0,0), (3,2,5,0,1,0), weights mod 6)
# dims from hess6_p1000000007.receipt.json; relaxation dims from the per-sector Singular runs (sectors.py and the hess6-quadrics cross-check);
# element certificates = verified exact Farkas certificates in certs/d6_p1000000007_k<key>_*.json (complement basis elements and random family points)

| sector | dim S_6 | dim K | dim J | non-J | relaxation affine dim | Farkas certs | verdict |
|---|---|---|---|---|---|---|---|
| (0, 0, 0) | 12 | 6 | 6 | 0 | 3 | 2 | trivial character: relaxation surface V0; excluded by reduced-system polynomial Farkas (Euler point D=6F inconsistent by hand) |
| (0, 0, 1) | 7 | 2 | 1 | 1 | 0 | 3 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 0, 2) | 8 | 2 | 0 | 2 | 0 | 4 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 0, 3) | 6 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (0, 0, 4) | 8 | 2 | 0 | 2 | 0 | 4 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 0, 5) | 7 | 2 | 1 | 1 | 0 | 3 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 1, 0) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 1, 2) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 1, 3) | 6 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (0, 1, 4) | 6 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (0, 1, 5) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 2, 0) | 8 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (0, 2, 2) | 8 | 2 | 0 | 2 | 0 | 4 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 2, 3) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 2, 4) | 8 | 2 | 0 | 2 | 0 | 4 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 3, 1) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 3, 5) | 5 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 4, 0) | 8 | 2 | 0 | 2 | 0 | 4 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 4, 3) | 5 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 5, 0) | 7 | 2 | 1 | 1 | 0 | 3 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (0, 5, 1) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 5, 4) | 5 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (0, 5, 5) | 6 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (3, 0, 0) | 6 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (3, 0, 1) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 0, 2) | 7 | 3 | 1 | 2 | 0 | 4 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (3, 0, 3) | 9 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (3, 0, 4) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 0, 5) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 1, 0) | 6 | 2 | 0 | 2 | 0 | 4 | whole sector family excluded (relaxation variety empty over F_p-bar => over Q-bar) |
| (3, 1, 1) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 1, 2) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 1, 3) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 1, 4) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 2, 1) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 2, 2) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 2, 3) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 3, 3) | 6 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (3, 4, 3) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 4, 4) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 4, 5) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 5, 1) | 5 | 1 | 0 | 1 | - | 1 | single non-Jacobian element: Farkas-excluded |
| (3, 5, 2) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 5, 3) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 5, 4) | 7 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |
| (3, 5, 5) | 6 | 1 | 1 | 0 | - | 0 | Jacobian element z_jF_i: Farkas certificate in reverify (certs/<Zj*F_Zi>.json) |

Totals: dim K_6 = 62, dim J_6 = 36, non-Jacobian = 26; sectors with kernel: 46 of 72.
