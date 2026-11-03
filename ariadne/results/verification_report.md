# ARIADNE verification report
Connectome: stanford_hardi (70 regions), onset P39, calibrated x0=-2.12

| Check | Result | Detail |
|---|---|---|
| Healthy network stays quiet | PASS | 0 regions recruited |
| Untreated seizure spreads | PASS | 69/70 recruited |
| Zero-strength drug equals baseline | PASS | 69 vs 69 |
| Result robust to noise seed | PASS | recruited across 3 seeds: [69, 69, 69] |
| Result converged in time step (dt 0.1 vs 0.05) | PASS | dt0.1 -> 69, dt0.05 -> 69 |
| Dose-response monotone: carbamazepine | PASS | recruited at strengths [0.0, 0.5, 1.0, 2.0, 4.0]: [69, 69, 69, 1, 1] |
| Dose-response monotone: lamotrigine | PASS | recruited at strengths [0.0, 0.5, 1.0, 2.0, 4.0]: [69, 69, 1, 1, 1] |
| Dose-response monotone: valproate | PASS | recruited at strengths [0.0, 0.5, 1.0, 2.0, 4.0]: [69, 69, 69, 1, 1] |
| Dose-response monotone: pregabalin | PASS | recruited at strengths [0.0, 0.5, 1.0, 2.0, 4.0]: [69, 69, 69, 69, 1] |
| Dose-response monotone: cenobamate | PASS | recruited at strengths [0.0, 0.5, 1.0, 2.0, 4.0]: [69, 1, 1, 1, 1] |
| Spread depends on onset site (not identical everywhere) | PASS | recruited from 4 other onsets: [69, 1, 69, 69] |

**11/11 checks passed.**

This verifies the simulator is internally consistent. It does NOT show the model predicts patient drug response.