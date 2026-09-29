# Model equations and limits

Only the induced ANN/bridge/balancer subgraph is modeled. ANN cells are traced
cables in the spatial model and single compartments in the point model. All
other selected cells are single compartments in both. Other cell types are
excluded; this is not the complete aboral organ.

For each compartment, with voltage v relative to a resting reference:

    C_i dv_i/dt = -(C_i/tau) v_i
                  + sum_j g_ij (v_j - v_i)
                  + sum_(k -> i) sign * I_syn * tanh(v_k/V_scale)
                  + I_ext,i(t)

The axial sum uses adjacent traced compartments. The synaptic sum uses one
contribution per measured directed contact, preserving autapses and polyadic
pre-to-post contacts. Counts are anatomical, not measured functional weights.
The odd tanh is a signed deviation-response approximation, not a nonnegative
firing rate or conductance-based receptor model.

Lengths are in µm, time in ms, voltage in mV, capacitance in pF, conductance in
nS and current in pA. CATMAID physical coordinates are converted from nm to µm.
Each node receives half the length of every incident edge, denoted l_i:

    C_i = 0.01 * C_m * 2*pi*r*l_i
    g_ij = pi*r^2*1e5 / (R_a*L_ij)

C_m is in µF/cm² and R_a in Ω·cm. The geometry represents lateral cylinder area,
without end caps or separately reconstructed soma surfaces. Three colocated
parent-child pairs are merged before conductances are computed.

The point model sums the spatial model's capacitances for each cell. Total leak,
injected current and synaptic contacts are matched. There is no separate weight
normalization that could bias the comparison.

| Assumption | Default |
|---|---:|
| Uniform cable radius | 0.25 µm |
| Axial resistivity | 150 Ω·cm |
| Specific capacitance | 1 µF/cm² |
| Leak time constant | 20 ms |
| Bridge/balancer cell capacitance | 0.5 pF |
| Synaptic current per contact | 0.05 pA |
| Synaptic voltage scale | 5 mV |
| Contact signs | All positive |
| Pulse | 1 pA, 20–40 ms |
| Timestep / duration | 0.5 / 160 ms |

The pulse site is the contact treenode with smallest source x coordinate on the
central ANN (ties broken by ID). It is a reproducible synthetic probe, not a
known sensory input. Passive membrane and axial currents use backward Euler;
bounded synaptic currents use the previous time step. Integration starts at rest.

Feedback removal disables 44 bridge-to-ANN contacts, leaving others unchanged.
The cut removes the central-tree axial edge that most evenly partitions traced
cable length; capacitance and synapses on both sides remain. It models neither
cell death nor wound healing. A point cell cannot represent this local cut.

Nine tests check source hashes, edge direction/counts, zero-length merging,
matched capacitance and axial charge conservation, zero input, the passive RC
cell-mean limit, the high-conduction point-model limit, timestep refinement, and
perturbation implementation. Passing these checks is not biological validation.

The 18-scenario sweep varies radius (0.1, 0.25, 0.5 µm), uniform sign (-1,+1)
and current (0,0.05,0.1 pA/contact). Zero-current scenarios repeat across signs
as controls. This is a scenario analysis, not a posterior or confidence interval.
Detailed convergence testing currently covers the default scenario, not every
possible parameter setting.
