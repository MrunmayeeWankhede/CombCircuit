# Live solver and interface

The browser kernel implements the MODEL.md equations at full source resolution.
Its passive matrix is a forest Laplacian plus diagonal membrane terms. Leaf-to-root
elimination computes the factorization, and root-to-leaf substitution solves each
implicit step in linear time. Synaptic currents use the preceding voltage and
remain directed; the synaptic circuit itself need not be a tree. The independently
integrated point circuit sums matched capacitances for each cell.

Changing parameters recomputes membrane values and the passive factorization;
voltage is retained. At a capacitance change, this means total stored charge can
change as part of the artificial intervention. A reset establishes a new clean
initial condition. One selected internal edge can be disabled, including any
edge selected in the viewer. Its conductance is set to zero; no compartment or
synaptic endpoint is deleted. The point circuit cannot express that internal cut.

The worker integrates fixed 0.5 ms steps. Rendering and requested playback speed
do not alter the timestep. Every retained trace sample is a calculated integration
step. Browser pacing may slow under load; it never jumps forward over equations.
The simulation pauses when the document is hidden.

Up to 12 probes can be selected from all 158 cells, all exported ANN compartments,
or three cell-type group means. Every original ANN treenode ID maps to its merged
compartment; IDs for colocated zero-length pairs share a compartment. Picking uses
the nearest projected segment within 12 screen pixels and chooses its nearest
endpoint, so overlapping paths should be separated by rotating or filtering.

A second model is integrated on the same inputs with 158 point compartments.
Point readouts for a local ANN probe read the single voltage of its owner cell.
Group and cell averages use capacitance weighting. No anatomy is generated for
the bridge or balancer cells, which are selected through the cell menu.

Exports retain up to 10,000 samples. An event log records reset settings, pulse
requests, parameter changes, probe sets, input locations and cuts. It is not a
serialized full-state checkpoint and there is no session-import feature. Reset
and probe changes clear the trace window as described in the interface.

Reference checks compare the entire JavaScript state to a SciPy sparse solve in
nine scenarios, including two cuts, negative input, individual bridge stimulation,
and parameter changes at 40 ms. These are implementation checks; they do not
validate electrical assumptions against the animal. Extreme parameter settings,
long-term behavior, unknown cable radii and missing syncytial loops remain open
scientific limitations.
