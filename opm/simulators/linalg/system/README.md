# System CPRW: choice of defaults

The CPRW pressure stage of `--linear-solver=system_cprw` has a few options whose
defaults (set in `setupSystemCPR`) were picked from runs on full Norne, with standard
wells and with every well converted to one segment per connection. The numbers are
total linear iterations for the whole run.

| Option | Values compared (one segment per connection) | Default |
|---|---|---|
| `well_weight_type` | cellavg 2569, cellblockavg 2604 | `cellavg` |
| `well_transfer` | classic 2558, no_prolongation 2569, full 2576 | `classic` |
| `well_coarse_diagonal` | contract_d 2569, row_sum 2580 | `contract_d` |

- `cellavg` and `cellblockavg` coincide for standard wells. `quasiimpes` behaves very
  badly for multisegment wells and should not become the default without re-checking
  them.
- The three `well_transfer` modes are within one iteration of each other for standard
  wells. The margin for `classic` is thin, and it was measured with an exact well
  solve, which nearly annihilates the well residual; restricting that residual may pay
  once the well solve is inexact.
- `contract_d` also brings the classic CPRW path itself from 2716 to 2646 on the same
  case.
- Summing a well's coarse column over every one of its blocks, rather than the top
  block alone, matters for multisegment wells: on Norne with one segment per
  connection, taking the top block alone loses every segment but the first -- most of
  the well -- making the coarse system far weaker than the classic cprw one.
