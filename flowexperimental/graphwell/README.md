# GraphWell — an entity-split multisegment well

GraphWell is an alternative formulation of the OPM multisegment well. It is
equivalent to the production `MultisegmentWell` for trees (matches it to
numerics), but is built so that **the flux unknowns are a separate entity from
the pressure/composition unknowns**. That separation is what gives it three
things the classic formulation cannot have cleanly:

1. a **transparent, face-based assembly loop**,
2. **consistent fixed-size automatic differentiation (AD)**, and
3. the structural freedom for **loops** (cyclic well topologies), because the
   number of flux DOFs is not tied to the number of segments.

This document explains the entity split, why it allows loops, and — in detail —
how the AD derivative slots are laid out.

---

## 1. Entities and degrees of freedom

There are two kinds of unknowns, living on two kinds of entities:

| entity | count | DOFs each | equation it carries |
|--------|-------|-----------|---------------------|
| **segment** | `nseg` | `NP` : pressure + (`NP−1`) phase fractions | `NP` component **mass-conservation** equations |
| **connection** ("face") | `nconn` | `1` : total volumetric flux | one **momentum / pressure-drop** equation (or, on the surface connection, the **well control** equation) |

`NP` = number of phases. A connection joins two segments (`down`, the deeper
end, and `up`, nearer the surface); one special **surface connection** sits on
the top segment and carries the control equation.

A **standard well** is just `1` segment + `1` surface connection. A **tree
multisegment well** has `nconn = nseg` (each segment has one outlet connection,
plus the surface connection replaces the top segment's outlet). A **loop** adds
extra connections **without** adding segments, so `nconn > nseg`.

---

## 2. Why the flux count is independent of the pressure/composition count

The whole well system is sized from `nseg` and `nconn` **independently**. The
total number of scalar unknowns is

```
flatSize = NP * nseg   +   nconn
           └ segments ┘     └ connections ┘
```

and the well Jacobian is a 2×2 `Dune::MultiTypeBlockMatrix`:

| block | shape | couples |
|-------|-------|---------|
| `Dss` | `nseg × nseg`  (blocks `NP×NP`) | segment mass eqns ↔ segment (p, fractions) |
| `Dsc` | `nseg × nconn` (blocks `NP×1`)  | segment mass eqns ↔ connection flux |
| `Dcs` | `nconn × nseg` (blocks `1×NP`)  | connection eqn ↔ segment (p, fractions) |
| `Dcc` | `nconn × nconn`(blocks `1×1`)   | connection eqn ↔ connection flux |

Nothing assumes `nconn == nseg`. The `Dcc` sparsity even couples **all
connections incident to a shared segment** (an `incident × incident` pattern),
so a junction where 3+ fluxes meet — the defining feature of a loop — is part of
the matrix structure, not a special case.

> The topology is stored as plain adjacency lists with **no single-outlet / tree
> assumption**, and upwinding is decided **per connection from the sign of that
> connection's own flux DOF** (`upwindSegment(c, Q)`), so the assembly never needs
> a global tree ordering. A unit test (`LoopBuildsAndFactorises`) adds a synthetic
> extra connection and checks the system still builds and factorises.

The only place a tree is still assumed is the **compatibility shim** that maps
GraphWell DOFs onto the production `MultisegmentWell` (so it can reuse MSW's
primary-variable / well-state / potential machinery). The core — topology,
equations, assembler — is loop-general; running a genuine loop only needs that
shim replaced by native primary-variable/well-state handling.

---

## 3. The AD layout in detail

All quantities are exposed as **fixed-size AD values** (`Opm::DenseAd::Evaluation`).
A single Eval type is used everywhere in the well assembly, with a derivative
layout chosen so that **any face term — which couples two segments and one flux —
fits in one Eval**. This is the key simplification over the production MSW, which
has to do a second "extra derivatives" evaluation for reverse flow because its
block can't hold both upwind and downwind derivatives at once.

### 3.1 The main Eval (`GEval`)

```
GEval = DenseAd::Evaluation<Scalar, NumResEq + 2*NP + 1>
```

with derivative slots:

```
slot index            meaning
─────────────────────────────────────────────────────────────
[0 ,  NumResEq)       reservoir cell variables (filled by perforation terms)
[RES, RES+NP)         "self"  segment : pressure at RES+0, fractions at RES+1..RES+NP-1
[OTH, OTH+NP)         "other" segment : pressure at OTH+0, fractions at OTH+1..OTH+NP-1
[QIDX]                connection flux
─────────────────────────────────────────────────────────────
RES  = NumResEq
OTH  = NumResEq + NP
QIDX = NumResEq + 2*NP
numDeriv = NumResEq + 2*NP + 1
```

For **3-phase black oil** (`NumResEq = 3`, `NP = 3`) that is **10 slots**:

```
 0  1  2 | 3  4  5 | 6  7  8 | 9
 └ res ─┘ └ self ─┘ └ other┘ flux
 R Rw Rg   p Fw Fg   p Fw Fg   Q
```

(`res` = the reservoir cell's conservation variables, used to build the
well→reservoir coupling; `self`/`other` = the two segments a face touches;
`flux` = the connection's own total-flux DOF.)

Mapping of the **main primary variables** to slots:

| primary variable | lives on | self-block slot | other-block slot |
|------------------|----------|-----------------|------------------|
| segment pressure `p`         | segment | `RES + 0`       | `OTH + 0` |
| segment fraction `F_k` (k=1..NP-1) | segment | `RES + k` | `OTH + k` |
| connection total flux `Q`    | connection | `QIDX`       | `QIDX` |

The "self" vs "other" distinction is purely positional: when a face term is
written for the connection's **down** segment, that segment's derivatives go in
the `self` block and the `up` segment's in the `other` block; a helper
(`shiftToOther`) relocates a value's derivatives from the self block to the other
block when the upwind segment is the `up` node. Because both blocks plus the flux
slot coexist in one Eval, **upwind and downwind derivatives are available
simultaneously** — no second evaluation.

The reservoir slots `[0, NumResEq)` are only non-zero for **perforation** terms,
which depend on the perforated cell's state (pressure, saturations, …). They feed
the `B`/`C` reservoir-coupling matrices.

#### The layout is a *superset* — each quantity fills only the slots it depends on

This is the key point to avoid confusion: the slot layout is the **union** of
everything any term might need. A given quantity is a *sparse occupant* of it.

- A **segment-local quantity** (segment pressure, a fraction, density, viscosity,
  the accumulation term) depends only on its own segment, so **only the `self`
  block is non-zero**; the `other` block and the `flux` slot are exactly zero.
  For example a segment density `rho_s` carries `d(rho_s)/d(p_s)` and
  `d(rho_s)/d(F_s)` in `[RES, RES+NP)` and nothing else.
- The **flux slot `[QIDX]`** is filled by exactly one primary variable — the
  connection flux `Q` (`connRate(c)` sets `QIDX = 1`, all other slots zero).
  Segment quantities never touch it; it only enters a segment equation through a
  *product*, e.g. the mass flux `q = Q · F(upwind)`.
- The **`other` block** is non-zero only in a genuinely two-sided term. In
  particular the **upwinded mass flux is one-sided**: `q = Q · F(upwind)` depends
  only on the *upwind* segment (placed in `self`) and on `Q` (`flux`), so its
  scatter is called with `other_col = -1` and the `other` block is unused. The
  `other` block earns its keep in the **momentum equation**, where
  `p_down − p_up` couples both segments.

`Role::Self` / `Role::Other` is therefore just a *placement choice*: it decides
whether a segment quantity's (self-block) derivatives are used as-is, or moved
into the `other` block by `shiftToOther`, so that a two-sided face term can hold
the `down` segment in `self` and the `up` segment in `other` at the same time. It
is the same quantity with the same derivatives, only repositioned.

Which slots each term actually populates:

| term (segment mass eqn unless noted) | self | other | flux | reservoir |
|--------------------------------------|:----:|:-----:|:----:|:---------:|
| accumulation                         | ✓ (this seg) | – | – | – |
| upwinded mass flux `q = Q·F(upwind)` | ✓ (= upwind seg) | – | ✓ | – |
| perforation inflow `cq_s`            | ✓ (this seg) | – | – | ✓ (cell) |
| momentum / pressure-drop (conn eqn)  | ✓ (down) | ✓ (up) | ✓ | – |
| control (surface conn)               | ✓ (top) | – | ✓ | – |

### 3.2 The compact control Eval (`CEval`)

The well control equation (BHP / surface-rate / THP) only depends on the **top
segment** and the **surface flux**, so it uses a smaller Eval:

```
CEval = DenseAd::Evaluation<Scalar, max(NP+1, 3)>
```

```
slot index   meaning
──────────────────────────────────────────
[0, NP)      top-segment DOFs : pressure at 0, fractions 1..NP-1
[QSlot]      surface-connection flux        (QSlot = NP)
──────────────────────────────────────────
```

(The size is padded to at least 3 so the shared `WellAssemble` /
`WellBhpThpCalculator` helpers — which are instantiated from size 3 — also cover
single-phase wells.)

### 3.3 How derivatives scatter into the matrix

Two generic helpers turn an assembled Eval into matrix entries — there is **no
manual offset arithmetic** anywhere in the term code:

`scatterMass(row_seg, comp, self_col, other_col, conn_col, term)` (a segment mass equation):

```
resSeg[row_seg][comp]                 +=  term.value()
Dss[row_seg][self_col ][comp][k]      +=  term.derivative(RES + k)   k = 0..NP-1
Dss[row_seg][other_col][comp][k]      +=  term.derivative(OTH + k)   k = 0..NP-1
Dsc[row_seg][conn_col ][comp][0]      +=  term.derivative(QIDX)
```

`scatterConn(conn_row, self_col, other_col, term)` (a connection momentum/control equation):

```
resConn[conn_row][0]                  +=  term.value()
Dcs[conn_row][self_col ][0][k]        +=  term.derivative(RES + k)   k = 0..NP-1
Dcs[conn_row][other_col][0][k]        +=  term.derivative(OTH + k)   k = 0..NP-1
Dcc[conn_row][conn_row ][0][0]        +=  term.derivative(QIDX)
```

Perforation terms additionally use the reservoir slots to fill the
reservoir-coupling blocks:

```
B[seg][perf][comp][k]  +=  cq_s[comp].derivative(k)          k = 0..NumResEq-1   (reservoir vars)
C[seg][perf][k][comp]  -=  cq_s[comp].derivative(RES + k)    k = 0..NP-1         (segment vars)
```

So the slot → matrix-block mapping is, in one glance:

| Eval slots | matrix block |
|------------|--------------|
| `[0, NumResEq)`      | `B` (reservoir → well), and `C` for the well → reservoir transpose |
| `[RES, RES+NP)`      | `Dss` / `Dcs` **self** block, and `C` |
| `[OTH, OTH+NP)`      | `Dss` / `Dcs` **other** block |
| `[QIDX]`             | `Dsc` / `Dcc` |

---

## 4. The assembly loop

The assembler iterates **entities**, not a tree:

- **over segments** — accumulation `(volume/dt)·(Sᵥ·frac − old)` and perforation
  inflow `cq_s` (drawdown × mobility, with crossflow handling), scaled by the
  well efficiency factor;
- **over connections** — one component flux `q = Q · F(upwind) · η` scattered into
  both adjacent segments with opposite sign (`+q` into `down`, `−q` into `up`),
  plus the connection's own equation:
  - internal connection → momentum: `p_down − p_up − Δp_hydro − Δp_friction/device − Δp_accel`,
    where the device term is the WSEGVALV / WSEGSICD / WSEGAICD drop when present;
  - surface connection → the well control equation.

Upwinding is frozen per Newton iteration from the sign of each connection's flux
DOF, which is what makes the loop case well-posed.

---

## 5. Linear system and reservoir coupling

The 2×2 multitype matrix is **flattened once** to a scalar BCRS matrix of size
`NP·nseg + nconn` and solved with UMFPACK for the Schur complement
(`flatSeg(s,eq) = NP·s + eq`, `flatConn(c) = NP·nseg + c`). The reservoir
coupling uses:

- `apply(x, Ax)` / `apply(r)` — the Schur operations `Ax −= Cᵀ·invD·B·x`, etc.;
- **CPRW** — `addWellPressureEquations` collapses the well to one well-pressure
  unknown coupled to the perforated cells (the segment-pressure row of `C`, and a
  perforation-averaged weight applied to `B`'s reservoir-pressure column).

---

## 6. Physics implemented (black oil)

Validated to numerics against the production `MultisegmentWell`:

- hydrostatic + friction + acceleration momentum (honours the `WELSEGS` `H__`/`HF-`/`HFA` flow model);
- BHP / surface-rate / **THP** control, including control switching;
- well efficiency factor (`WEFAC`);
- flow-control devices **WSEGVALV**, **WSEGSICD**, **WSEGAICD** (per-connection device drops + `SPRD*` reporting);
- implicit IPR;
- **CPRW** well-pressure coupling.

---

## 7. Status and the loop outlook

The GraphWell currently derives from `MultisegmentWell` and reuses its
primary-variable / well-state / potential / THP-operability / segment-setup
infrastructure; that reuse layer is the only part that assumes a tree (via the
`Q_c = −WQTotal_{down(c)}` bijection). The equation core — entity split,
independent flux DOFs, face-based AD assembly, and the `Dcc` shared-segment
coupling — is already loop-general. Realising an actual looped well needs that
reuse layer replaced with native primary-variable/well-state handling (the
deferred "derive directly from `WellInterface`" step); MPI-distributed wells are
also deferred.
