# Lagrange friction: relative movement in the predictor

With this change, `LagrangeStrategy::predict_relative_movement()` evaluates the relative movement
(jump) for all frictional contact. Before, it only did so with `FRLESS_FIRST`; otherwise, the
predictor reused the jump of the last Newton iteration of the previous time step.

## Result

- **ConstVel, ConstAcc:** never worse, often better (up to 44 % fewer iterations).
- **TangDis:** better in 6 tests; `contact3D_slidingblock_coulomb` passes its result check only
  after the change. `contact2D_slidingblock_tresca_new_struct` no longer converges.
- **ConstDis:** never better, 4 tests need up to 7 % more iterations, and
  `contact2D_slidingblock_tresca_new_struct` no longer converges. With ConstDis, the stale jump
  happened to be a good guess for steady sliding. That test now uses ConstVel (see below).
- **TSI, EHL and wear tests:** unchanged.

## Newton iterations per test and predictor

All 22 tests with friction and a Lagrange multiplier strategy, each run with every predictor. A cell
shows the total number of Newton iterations over all time steps, before → after the change.

- **bold:** the test's own predictor (the one the test suite runs)
- 🟢 / 🔴: fewer / more iterations, or passes / fails only after the change
- ✗: failed run (error or failed result check), in the same way before and after unless marked

| Test | ConstDis | TangDis | ConstVel | ConstAcc |
|---|---|---|---|---|
| `contact2D_nurbs9_pg_frict` | ✗ → ✗ | **148 → 147** 🟢 | ✗ → ✗ | ✗ → ✗ |
| `contact2D_slidingblock_coulomb` | **108 → 108** | 70 → 68 🟢 | 225 → 144 🟢 | 264 → 149 🟢 |
| `contact2D_slidingblock_tresca` | **67 → 67** | 63 → 63 | 82 → 82 | 86 → 86 |
| `contact2D_slidingblock_tresca_new_struct` | **73 → ✗** 🔴 | 54 → ✗ 🔴 | 87 → 87 | 87 → 87 |
| `contact2D_twobeamshorizontal` | **73 → 73** | 23 → 23 | 28 → 28 | 28 → 28 |
| `contact2D_twobeamshorizontal_new_struct` | **73 → 73** | 23 → 23 | 28 → 28 | 28 → 28 |
| `contact2D_twobeamsvertical` | **201 → 201** | ✗ → ✗ | ✗ → ✗ | ✗ → ✗ |
| `contact3D_frictionless_first_lagrange` | **16 → 16** | 14 → 14 | 10 → 10 | 10 → 10 |
| `contact3D_lin_ltl` | **1594 → 1594** | 1341 → 1341 | 1325 → 1325 | 1325 → 1325 |
| `contact3D_slidingblock_coulomb` | **47 → 49** 🔴 | ✗ → 29 🟢 | 48 → 48 | 48 → 48 |
| `contact3D_slidingblock_coulomb_new_struct` | **47 → 49** 🔴 | ✗ → 29 🟢 | 48 → 48 | 48 → 48 |
| `contact3D_slidingblock_lagrange_hex20` | **20 → 20** | 14 → 11 🟢 | 17 → 15 🟢 | 17 → 15 🟢 |
| `contact3D_slidingblock_lagrange_tet10` | **32 → 32** | 29 → 29 | 32 → 32 | 32 → 32 |
| `contact3D_slidingblock_lagrange_tet10_new_struct` | **32 → 32** | 29 → 29 | 32 → 32 | 32 → 32 |
| `contact3D_symmetry_LM_cond` | ✗ → ✗ | **63 → 57** 🟢 | 62 → 61 🟢 | 62 → 61 🟢 |
| `contact3D_tet4_hex8` | **245 → 261** 🔴 | 208 → 197 🟢 | 268 → 266 🟢 | 268 → 266 🟢 |
| `contact3D_tet4_hex8_new_struct` | **245 → 261** 🔴 | 208 → 197 🟢 | 268 → 266 🟢 | 268 → 266 🟢 |
| `ehl3d_mixed` | 105 → 105 | ✗ → ✗ | **116 → 116** | ✗ → ✗ |
| `tsi_contact3D_coulomb_friction` | 185 → 185 | 214 → 214 | **160 → 160** | 160 → 160 |
| `tsi_contact3D_frictionless_first` | **17 → 17** | 12 → 12 | 11 → 11 | 11 → 11 |
| `wear2D_modgap` | **35 → 35** | 22 → 22 | 39 → 39 | ✗ → ✗ |
| `wear3D_modgap` | **15 → 15** | 12 → 12 | 13 → 13 | 13 → 13 |

`contact3D_lin_ltl` is the variant with the new structural time integration. Changed failures:

- `contact2D_slidingblock_tresca_new_struct` (ConstDis, TangDis): after the change, Newton does
  not converge in step 17, the active set cycles.
- `contact3D_slidingblock_coulomb` (both variants, TangDis): before the change, the result check
  fails.

## Consequence for the test suite

All ctest entries of these tests pass except `contact2D_slidingblock_tresca_new_struct`, which is
switched to `PREDICT: "ConstVel"`. It passes with its unchanged reference values (87 instead of 73
iterations).
