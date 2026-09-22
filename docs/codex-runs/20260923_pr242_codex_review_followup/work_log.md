# PR #242 Codex review follow-up: Issue #61 raw reanalysis

## Immutable input and local synchronization

The source campaign was not changed:

`/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-09-06/220250/parallel/issue61-2015-1tau__3654804140d9/campaign`

The three raw `step_summary.csv` inputs were synchronized to the corresponding local campaign conditions and SHA-256 verified:

| condition | SHA-256 |
| --- | --- |
| `project_torque_1em21` | `ece66917fef20475c3020627d418862aff95ee83ff148c2cafe4704f2a9ecd24` |
| `project_torque_2p5em20` | `f983a40d7907cb20ff73ba8162bd72dbe9d86aff60b8a820c98751b414baf6df` |
| `project_torque_1em19` | `e9edcc7e9b6c92b2b84475f0b8be685b2b3ee34eb9a054fea60ebd1cf98eca82` |

## Corrected result

The updated streaming analysis was run on cs10 without a simulation restart. It wrote
`outputs/2026-09-23/023008/analysis/issue61_2015_1tau_corrected/`; the two result files were synchronized locally and their SHA-256 values matched cs10.

| file | SHA-256 |
| --- | --- |
| `issue61_summary.csv` | `528d28bf2f744656e225cb51600935f7dcf4bc8835dccbf6aaa9e40b17ff74e9` |
| `issue61_decision.json` | `a1269b96bfbe5136f21b98694d13e42f8a962e4fd941d06349da61fe28b5537b` |

All three conditions remain strict FAIL. The first observed crossing is
`max_motor_torque_balance_residual_ratio` at step 0: `1e-5 s` (`1e-21 N m`),
`4e-7 s` (`2.5e-20 N m`), and `1e-7 s` (`1e-19 N m`). Helix-pitch violations also occur, but later; the locked threshold contract was not changed.

Reservation 8 (`8278c82`) remains running and reservation 6 remains queued; neither was started, stopped, or altered by this reanalysis.
