# First results

This is a computational implementation milestone, not a biological discovery.

The default model injects 1 pA into a reproducibly selected location on ANN Q1-4 from 20 to 40 ms. It assumes a uniform 0.25 µm radius and excitatory synaptic contacts. Both representations have matched total capacitance and input current.

| Readout | Equipotential peak (mV) | Spatial peak (mV) |
|---|---:|---:|
| ANN Q1-4 | 0.359292 | 0.357629 |
| ANN Q1Q2 | 0.010176 | 0.009841 |
| ANN Q3Q4 | 0.023360 | 0.020830 |
| bridge | 0.123596 | 0.113291 |
| balancer | 0.032340 | 0.033048 |
| stimulated compartment | 0.359292 | 0.886557 |

The local response differs more than the central cell average in this illustrative scenario. This motivates a measurement-design question: could local readouts distinguish plausible mechanisms when averages cannot? It does not establish what the real animal does.

The balancer values are modeled electrical drive, not ciliary beat frequency. Neither model has been fitted to functional measurements. A larger spatial response is not evidence of greater biological accuracy.

The 18-scenario sensitivity sweep varies radius, all-contact sign and current amplitude. It is a scenario analysis, not confidence bounds. With synaptic currents off, the mean responses agree to numerical precision. With them on, differences depend strongly on assumptions.

Nine numerical tests passed, including source hashes, contact direction, charge conservation, the passive RC limit, the high-conduction limit and timestep refinement. Those checks test implementation rather than biological validity.

The next step is to resolve geometry and physiological constraints, then find a held-out perturbation that separates plausible calibrated models. See STUDY_PLAN.md.