# 3D mortar: ambiguous Delaunay triangulation of clip polygons

## Summary

`Mortar::Coupling3d::delaunay_triangulation()` splits a clip polygon into integration triangles.
If the vertices of the polygon lie on one circle, e.g. for a rectangle or a symmetric trapezoid,
both diagonals are valid Delaunay triangulations, and the circumcircle check, which has no
tolerance, decides between them. The integrands are not polynomial in general, so the two
diagonals give integrals that differ by the quadrature error. In contact, the triangulation is
redone in every Newton iteration, and in symmetric problems the diagonal can flip in every
iteration. Newton then alternates between two mirror-image states and does not converge.

**Status:** parked. The fix on this branch (fall back to the center-based triangulation near
cocircular configurations) works for the reproducer, but depends on a tolerance and does not
remove the underlying discontinuity, see [Alternatives](#alternatives). **Workaround:** set
`TRIANGULATION: "Center"` in `MORTAR COUPLING` if Newton cycles in a symmetric contact problem.

## Where the clip polygon is triangulated

`Mortar::Coupling3d::evaluate_coupling()` builds the clip polygon of a slave and a master element
(see `fix-mortar-clipping-hull` for that part) and then calls `triangulation()`:

1. `TRIANGULATION: "Delaunay"` (the default): `delaunay_triangulation()`. If it fails (no
   well-defined circumcircles), the center-based triangulation is used as a backup.
2. `TRIANGULATION: "Center"`: `center_triangulation()` connects the area centroid of the polygon
   with all its edges, so a polygon with N vertices gives N triangles instead of N−2. A triangle
   is used directly as one cell.

`delaunay_triangulation()` builds the triangulation incrementally: it starts with the triangle of
the first three vertices, adds the remaining vertices one by one, and checks for every triangle
whether another vertex lies inside its circumcircle (`dist < radius1`, without tolerance). The
tolerance `tol` is only used afterwards, to make close good/bad decisions of different triangles
consistent with each other.

The Gauss points of the cells are projected onto the slave and master elements along the constant
normal of the auxiliary plane. The gap uses the interpolated nodal normals, normalized to length 1
(`CONTACT::Integrator::gap_3d()`).

The volume mortar coupling (`Coupling::VolMortar`) has its own `delaunay_triangulation()` and is
not considered here.

## When Newton cycles

For a quadrilateral clip polygon, the code evaluates one of two residuals: R_A(u) if the
circumcircle check selects diagonal A, R_B(u) if it selects B. The switch happens on the set S of
configurations in which the polygon is exactly cocircular, and there the residual jumps by the
difference of the quadrature errors of the two triangulations. Each branch has its own root u_A*
and u_B*, both close to the exact solution:

- If each root lies on its own side of S, Newton picks one side and converges.
- If each root lies on the other side of S, R has no root at all. Newton goes towards u_A*, lands
  on the B side, goes towards u_B*, lands on the A side, and so on.

All four of the following conditions are needed for Newton to fail:

1. **A discontinuity in the residual**, here the switch of the diagonal. The two triangulations
   must give different integrals, which needs non-polynomial integrands, e.g. distorted faces
   (from a different lateral expansion of the two bodies), curved contact surfaces (via the
   normalized normal in the gap), or a contact or stick/slip boundary inside a polygon (kink in
   the integrand). With equal materials, or for aligned identical cubes, Newton converges although
   the decision flips.
2. **Something that pins the solution onto the discontinuity.** Normally the solution crosses S
   at one load level and leaves it again. A symmetry keeps it there: if a slave/master element
   pair is cut into mirror halves by a symmetry plane, its clip polygon stays a symmetric
   trapezoid, which is always cocircular. With an odd number of elements, this is the middle
   pair. With an even number, the symmetry plane lies on element boundaries and nothing happens.
3. **The wrong sign.** The mirror symmetry maps diagonal A to B, so the two discrete roots lie at
   ±δ around S, where δ is the asymmetry caused by the quadrature. Whether each root lies on its
   own side or on the other side is a property of the problem. In the reproducer it is the other
   side: the asymmetry left by one diagonal selects the other diagonal in the next iteration.
4. **A jump above the tolerances.** Newton alternates between two states whose residuals differ
   by the jump, so the residual floor is bounded by the quadrature error. It grows with the
   distortion of the faces and decreases with mesh refinement for smooth integrands, but only
   slowly if the integrand has a kink. Newton fails only if the tolerances are below this floor.

Friction, sliding, the predictor and the step size do not matter.

### Pseudo-2D models in plane strain

With one element through the thickness and the displacements in the thickness direction fixed,
all contact faces are flat rectangles (a 2D segment extruded in thickness direction). Every clip
polygon is then a rectangle, i.e. cocircular, and the mirror symmetry through the mid-plane pins
every slave/master pair onto S, not only a middle one. The projection onto the flat faces is
affine, so the shape function products in D and M are integrated exactly by both diagonals. But
the gap uses the normalized interpolated normal:

- **Flat contact** (all nodal normals equal): all integrands are polynomial, so there is no jump.
  It cannot happen.
- **Curved contact** (e.g. a cylinder on a block): the normalization makes the gap rational in
  the in-plane direction, and the two diagonals give front/back mirror-image results that differ.
  The jump scales with the change of the normal across an element and should be small for
  resolved curvatures.
- **Contact or stick/slip boundaries** inside a polygon: kink in the integrand, larger jump.

Not tested: the plane-strain variant of the reproducer (expected to converge) and a coarse
cylinder on a block (expected to flip in every pair, with a nonzero residual floor).

## Consequences

### Ambiguous triangulations in existing tests

`delaunay_ambiguous_diagnostic.patch` (apply on main) prints
`AMBIGUOUS DELAUNAY TRIANGULATION: ...` whenever a clip polygon vertex lies within a relative
distance of 1e-6 of a circumcircle (the tolerance of the fix below), together with the element
pair and the side the vertex was put on. Comparing these lines between Newton iterations shows
whether the triangulation flips. On main, for the 194 inputs with 3D mortar coupling that do not
set `TRIANGULATION: "Center"` (without beam, embedded mesh, FBI and volume mortar tests), it
reports:

| Inputs | With reports | Same vertex on both sides during the run | Reports |
|---|---|---|---|
| Contact (`contact3D_*`, `roughcontact3d_*`, `wear3D_modgap`, `tsi_contact3D_conduction`, `fsi_pw_mono_fs_ga_ga_contact`, `fsci_simp_salz_syst`) | 49 | 41 | 56345 |
| Meshtying, FSI mortar, S2I mortar (`meshtying3D_*`, `f3_*_meshtying_*`, `fsi_*_mtr*`, `*_nurbs`, `sti_*_mortar_*`) | 54 | 8 | 36374 |

Only the reproducer below fails. The reports fall into two groups:

- **Exactly cocircular** (relative distance up to 1e-13): rectangles from matching or aligned
  structured meshes, where round-off decides. This is all of meshtying and the simple contact
  patch tests. The residuals of e.g. `contact3D_lin_*`, `contact3D_roundrobin_ghost` and
  `contact3D_patch_linstatic` with Delaunay and with the fallback differ only at about 1e-13, so
  both diagonals give practically the same integrals.
- **Nearly cocircular** (1e-9 to 1e-6): deforming or sliding problems (`contact3D_slidingblock_*`,
  the hex27/tet10/hex20 penalty tests, `contact3D_sp_*`, `fsi_pw_mono_fs_ga_ga_contact`,
  `fsci_simp_salz_syst`). These are not ambiguous: the circumcircle check is well-defined there.
  They are only reported because the band of 1e-6 is wide.

None of these tests has a symmetry plane through an element pair together with distorted faces
and tight tolerances, so condition 2, 3 or 4 is missing.

### Newton iterations cycle and do not converge

The reproducer `contact3D_delaunay_symmetric_two_hex8` has two frictionless HEX8 cubes. The
stiffer upper cube (slave) is offset by 0.3 in x and pressed into the lower one (master). The
problem is symmetric in z, and the faces in z are free, so the overlap is a symmetric trapezoid
whose shape depends on the different lateral expansion of the two cubes.

On main (2 ranks, as in the test), the diagnostic reports the same element pair in every Newton
iteration of every step, and the side of the vertex alternates in every iteration. In step 10 the
vertex is at a relative distance of ±1.6e-11 from the circumcircle, about 1e5 times the machine
precision. Each step ends at a residual floor that grows with the deformation, from 9.8e-12 in
step 1 to 4.8e-9 in step 10. In step 10 the increment norm stays at 1.14e-10 > `TOLDISP` 1e-10,
so the step does not converge:

```
it 0:   ||F|| = 3.76607e+00  dx = 0.00000e+00   vertex inside
it 1:   ||F|| = 4.36695e-02  dx = 2.48021e-02   vertex outside
it 2:   ||F|| = 1.49259e-05  dx = 4.60872e-04   vertex inside
it 3:   ||F|| = 4.84812e-09  dx = 4.08710e-07   vertex outside
it 4:   ||F|| = 4.84812e-09  dx = 1.13761e-10   vertex inside
it 5:   ||F|| = 4.84812e-09  dx = 1.13761e-10   vertex outside
...
it 19:  ||F|| = 4.84812e-09  dx = 1.13761e-10   vertex outside
it 20:  ||F|| = 4.84812e-09  dx = 1.13761e-10   vertex inside  (Failed!)

Nonlinear solver did not converge. Aborting since DIVERCONT is set to stop.
```

The residual and increment norms stay constant while every entry of the increment flips its sign
from one iteration to the next. With `TRIANGULATION: "Center"`, or with the upper cube shifted by
0.001 in z, all steps converge in 4 iterations.

The reproducer is mild: the floor is about 1e-9 relative to the initial residual, and with
`TOLDISP` 1e-8 it would pass. It uses one element per body, i.e. the coarsest possible mesh.
Coarse meshes, larger distortions and non-smooth integrands give larger floors. In
`contact3D_nitsche_hex8_sym_coulomb`, switching from Delaunay to the fallback changes the first
residual by 25 % (0.475 to 0.594), so the jump can be far larger than in the reproducer.

## Fix on this branch

If a vertex lies within `MORTARDELAUNAYTOL * radius` (1e-6) of a circumcircle,
`delaunay_triangulation()` returns `false`, so the existing center-based backup is used. For a
symmetric polygon, the centroid lies on the symmetry axis and the center-based triangles are
mirror images of each other, so this triangulation does not make the solution asymmetric.

This is only done for contact (`CONTACT::Coupling3d::center_triangulation_if_ambiguous()`), whose
mortar terms are evaluated in every Newton iteration. Meshtying evaluates them once, so it keeps
the Delaunay triangulation and its reference results.

The fix moves the discontinuity from S to the edge of the tolerance band, which no symmetry pins
the solution onto. It does not remove it, and the band has to be wider than the asymmetry δ of
nearly symmetric problems, which is not known in advance. The data would also allow a band of
1e-10 (exact ties up to 1e-13, the reproducer's cycling states at 1.6e-11, unambiguous polygons
from 1e-9), but this was not tried.

### Result changes

The fallback is active in all 48 other contact inputs from the table above. All of them still
pass with their reference results. The largest change of a result value, apart from round-off,
is about 4.5e-11 (`contact3D_slidingblock_coulomb`), and the number of Newton iterations is
unchanged, except for `contact3D_quad_tet10` (6 instead of 7 iterations in step 1). In the hex27
and hex20 penalty tests, the fallback is active (632 times in `contact3D_pwlin_penalty_hex27`),
but the results are bit-identical, presumably because the affected pairs only contribute to
inactive slave nodes (not verified). Meshtying, FSI mortar and S2I mortar are not affected.

The reproducer and all variants I tried (frictionless and Coulomb, pressing and sliding, TSI)
converge with the fix, with the same iteration counts as with Center: 4 iterations in every step,
down to a residual of about 1e-14. New test:

- `tests/input_files/contact3D_delaunay_symmetric_two_hex8.4C.yaml`: the reproducer above, which
  does not converge on main.

## Alternatives

### Why not make the Delaunay triangulation itself unambiguous

- **Exact predicates** (e.g. an exact in-circle test) decide exactly for the given floating-point
  coordinates. But the cycling iterates are not cocircular: they lie at ±δ from the circumcircle,
  where the check is well-defined. The decision still flips.
- **A tie-break at round-off level**, e.g. by global vertex IDs, has the same problem.
- **A fixed tie-break within a tolerance band**, e.g. always the diagonal through the vertex with
  the smallest global ID, only works if the band is wider than δ, and makes the solution
  asymmetric by δ. A tolerance in the circumcircle check alone was tried: Newton then alternated
  between two mirror-image states with an asymmetry just above the tolerance.

Any rule that chooses between two diagonals has a switching set somewhere. What matters is
whether the symmetries of the problem can pin the solution onto it.

### Options without a tolerance

| | Cycling on the diagonal | Tolerance | Continuous | Frame-invariant | Cells |
|---|---|---|---|---|---|
| Delaunay (main) | yes, for symmetric problems | — | no | yes | N−2 |
| Delaunay + fallback (this branch) | only if the band is narrower than δ | yes | no | yes | N−2 / N |
| Center for contact | no | no | yes, apart from the triangle shortcut | yes | N |
| Fan from the vertex that is extreme in a generic direction | not by coordinate-plane symmetries | no | no | no | N−2 |
| Freeze the diagonal within a time step | not within a step | no | within a step | yes | N−2 |

- **Center for contact:** the area centroid moves continuously with the vertices, so the
  triangulation is continuous, and it is symmetric for symmetric polygons. A new vertex appearing
  on an edge does not move the centroid either, but the shortcut for triangles (1 cell instead of
  3 around the centroid) makes the transition from 3 to 4 vertices discontinuous. It would have
  to be removed as well. Costs: about twice as many quadrature points, many changed reference
  results, a changed default.
- **Fan from an extreme vertex:** triangulate from the vertex with the largest g·x, with a generic
  direction like g = (1, √2, π). It only switches if an edge is perpendicular to g, which a mirror
  symmetry about a coordinate plane cannot pin. In the reproducer, the choice is decided by a wide
  margin. But the result depends on the orientation of the model, and the triangles can be thin.
- **Freeze within a time step:** decide the diagonal at the beginning of the step and keep it
  during the Newton iterations. Needs a cache of the choice per slave/master pair and vertex
  identities, since the coupling objects are rebuilt in every evaluation, and a rule for clip
  polygons that change their vertices during the step. The result becomes slightly
  history-dependent.

## Not changed yet: meshtying

Applying the fallback to meshtying as well changes 7 meshtying tests: by 1e-9 to 4e-7 for
Lagrange elements (`meshtying3D_structure_penalty_redist_static`,
`meshtying3D_structure_quad_const`, `meshtying3D_structure_quad_lin`, each also `_new_struct`
where present), but by 4e-5 (`thermo3D_meshtying_nurbs`) and 2 % (`meshtying3D_nurbs_dual`) for
NURBS. The NURBS quadrature seems to be poorly resolved there. For meshtying the triangulation is
only done once, so it cannot flip between iterations, but the result still depends on round-off.

## Open points

- Decide between the fallback and one of the options without a tolerance. Making `Center` the
  default for contact is the cleanest one (43 inputs already set it explicitly), but needs input
  from the mortar maintainers.
- Test the claims about plane strain and mesh refinement: plane-strain reproducer, coarse
  cylinder on a block, and the reproducer with n×n×n elements per cube (n = 1, 3, 5).
- A separate bug in `Mortar::sort_convex_hull_points` (dropped clip polygon corners for nearly
  collinear points) is fixed on branch `fix-mortar-clipping-hull`.
