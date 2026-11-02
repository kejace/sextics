events: {'coker': 21, 'e1space': 21, 'relax': 3, 'relax-reconstruct': 1, 'confirm': 1, 'msolve-export': 18, 'msolve': 11, 'k8-e1': 1, 'k8-fixedD-symmetric': 1, 'wave2-search-start': 1, 'k8-fixedD-full': 1, 'famrelax': 2, 'msolve-family-export': 2, 'msolve-family': 2}

### E1 spaces (exact mod p)
| d | sec | dim S | dim J | dim R | E1-admissible | non-J |
| 4 | 00 | 48 | 0 | 48 | 0 | 0 |
| 4 | 01 | 38 | 0 | 38 | 0 | 0 |
| 4 | 10 | 22 | 0 | 22 | 0 | 0 |
| 4 | 11 | 18 | 0 | 18 | 0 | 0 |
| 5 | 00 | 88 | 3 | 85 | 3 | 0 |
| 5 | 01 | 78 | 2 | 76 | 2 | 0 |
| 5 | 10 | 48 | 1 | 47 | 1 | 0 |
| 5 | 11 | 38 | 0 | 38 | 0 | 0 |
| 6 | 00 | 158 | 14 | 144 | 14 | 0 |
| 6 | 01 | 138 | 12 | 126 | 12 | 0 |
| 6 | 10 | 88 | 6 | 82 | 6 | 0 |
| 6 | 11 | 78 | 4 | 74 | 4 | 0 |
| 7 | 00 | 258 | 45 | 213 | 45 | 0 |
| 7 | 01 | 238 | 40 | 198 | 40 | 0 |
| 7 | 10 | 158 | 23 | 135 | 23 | 0 |
| 7 | 11 | 138 | 18 | 120 | 18 | 0 |
| 8 | 00 | 413 | 112 | 301 | 113 | 1 |
| 8 | 01 | 378 | 104 | 274 | 104 | 0 |
| 8 | 10 | 258 | 64 | 194 | 64 | 0 |
| 8 | 11 | 238 | 56 | 182 | 56 | 0 |

### Hessian relaxation, exact (msolve)
| d | sectors | family | K | polys | empty |
| 6 | 00 | J | 14 | 296 | False  |
| 5 | 00 | J | 3 | 23 | True  |
| 5 | 01 | J | 2 | 19 | True  |
| 5 | 10 | J | 1 | 8 | True  |
| 6 | 01 | J | 12 | 294 | True  |
| 6 | 10 | J | 6 | 234 | True  |
| 6 | 11 | J | 4 | 66 | True  |
| 6 | e | J | 18 | 418 | True  |
| 6 | o | J | 18 | 473 | True  |
| 7 | 11 | J | 18 | 426 | True  |
| 7 | 10 | J | 23 | 427 | True  |
| 6 | 00 | J | 14 | 296 | True  |
| 7 | 11 | J | 18 | 426 | True  |

### Relaxation LM landscapes
| d | sectors | family/affine | K | search dim | m | restarts | best | median |
| 5 | ['00', '11'] | None | 3 | None | 744 | 10 | 5.98e-01 | 5.99e-01 |
| 6 | ['00'] | None | 14 | None | 492 | 1 | 1.49e-15 | 1.49e-15 |
| 6 | ['00'] | None | 14 | 13 | 492 | 30 | 5.60e-01 | 6.09e-01 |
| 6 | ['00'] | J | 14 | 13 | 492 | 30 | 5.60e-01 | 6.09e-01 |
| 8 | ['00'] | JK | 113 | 103 | 274 | 40 | 2.75e-03 | 2.77e-03 |

### Koszul-reduced full-E2 varpro
| d | sectors | family | K | search dim | rows | n_a | restarts | best rel | median | time |

### Exact confirmations
- d=6 spurious (F | D) rank=None nullity=None note=relax d=6 secs=00 rel=1.49e-15 D=6*X0**6 + 30*X0**3*X1*Y0*Y1/7 + 30*X0*X1*Y0**3*Y1/7 + 30*X0*X1*Y0**2*Y
- k8-e1: {"e1": true, "info": {"rank": 5052, "rank_aug": 5052, "rows": 5793, "cols": 6006}, "seconds": 82.5}
- k8-fixedD-symmetric: {"status": "inconsistent", "seconds": 36.1, "prime": 2147483647, "rank": 3781, "nullity": 475, "rows": 22016, "cols": 4256, "full_cols": 7497, "subgroup_size": 4, "group_size": 4, "d": 8}
- k8-fixedD-full: {"seconds": 195.3, "consistent": false, "rank": 6703, "nullity": 794, "rows": 40301, "cols": 7497, "d": 8}
- relax-reconstruct: {"d": 6, "rel": 1.4934890565976107e-15, "exact_relaxation": true, "g": {"0:(1, 0, 0, 0, 0, 0)": "1", "1:(0, 1, 0, 0, 0, 0)": "1", "2:(0, 0, 1, 0, 0, 0)": "1", "3:(0, 0, 0, 1, 0, 0)": "1", "4:(0, 0, 0, 0, 1, 0)": "1", "5:(0, 0, 0, 0, 0, 1)": "1"}}
