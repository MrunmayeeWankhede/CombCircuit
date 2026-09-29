# CombCircuit Live 0.2 — video transcript

Synthetic narration; measured anatomy and assumed physiology are distinguished throughout.



## 1. A live circuit you can experiment with



Welcome to Comb Circuit Live, version zero point two. This walkthrough explains the live simulator, its controls, every adjustable parameter, and how to interpret your experiments.

The old explorer replayed saved curves. This version calculates new voltages as you run it, including responses to inputs and parameter changes you choose.

Unzip the project and open LIVE dot H T M L in a browser. No account, Python installation, or internet connection is needed to use the simulator.

The anatomy is measured, but the electrical properties are hypotheses. Live computation does not mean we have recorded activity from a living animal.



## 2. The circuit inside the comb jelly



The organism is the comb jelly Mnemiopsis leidyi. The selected circuit comes from its aboral organ, which is associated with balance and ciliary coordination.

The reconstruction contains three aboral nerve net neurons, abbreviated A N N, twenty six bridge cells, and one hundred twenty nine balancer cells.

That is one hundred fifty eight cells, including non-neuronal balancer cells, with two hundred eighty directed anatomical contacts. It is not the whole animal.

The original anatomical work belongs to Jokura and colleagues. Our contribution is an exploratory electrical modeling workbench built on their data.



## 3. Branches, cells and compartments



The branching picture shows only the three reconstructed A N N arbors. Q one through four spans all four anatomical quadrants; the other two arbors each span two quadrants.

These neurons have syncytial morphology, including continuous branched structures and multiple nuclei. The visible lines are traced skeletons, not thousands of separate neurons.

The spatial model divides the arbors into sixteen thousand four hundred forty six electrical compartments. Bridge and balancer cells each use one compartment, giving sixteen thousand six hundred one altogether.

Line crossings in a screen projection are not proof of a connection. The exported skeleton trees may also omit internal reconnecting loops, so cuts must be interpreted cautiously.



## 4. Two models, calculated together



The spatial model allows voltage to differ between locations within one A N N neuron. Current spreads through its reconstructed internal paths and leaks through its membrane.

The point model, also called equipotential, gives each cell just one voltage. It cannot represent a local gradient within that cell.

Both models receive the same injected current and anatomical contacts, and have matched total cell capacitance and leak. Solid graph lines show the spatial model; dashed lines show the point model.

A local probe can disagree strongly between models even when the whole-cell averages look similar. More anatomical detail alone does not establish which response is biologically correct.



## 5. Run, pause, step and reset



Press Run to begin, and press Pause to stop at the current state. Step advances ten milliseconds of simulated time without starting continuous playback.

With the demo option checked, a one picoampere pulse is scheduled from twenty to forty milliseconds after a reset. The pulse uses the currently chosen input site and pulse settings.

The speed selector controls requested simulated time per wall-clock second. One times biological time requests one thousand simulated milliseconds each second, but slower computers may achieve less.

The clock always shows simulated time. The solver does not skip numerical steps to catch up. Reset clears voltages, time, trace history and the event log, while keeping your current probes and settings.



## 6. Choose your own local probes



Set Click action to Add local probe, then click a branch. A colored ring marks the selected compartment, and a matching trace appears in the measurement panel.

The click selects a nearby traced compartment in the current screen projection. Rotate or zoom to separate overlapping branches before choosing a location.

You can display up to twelve arbitrarily chosen probes together. Remove one with its cross button to free a channel. The limit is for readable display, not six fixed biological measurements.

Adding or removing a probe starts a new trace window without resetting the electrical state. Newly added probes cannot show activity that was never recorded for them.



## 7. Read any cell, group or treenode



The Any cell menu includes all one hundred fifty eight cells. Choose an individual bridge or balancer cell and press Add cell mean to inspect its response.

For an A N N cell, the mean averages its compartments with weights proportional to capacitance. Larger membrane contributions count more than small pieces.

The group menu adds means across bridge cells, balancer cells, or all A N N neurons. A local node probe instead reads one specific compartment.

You can also enter an original A N N treenode identifier. In the point comparison, every local probe within a particular cell necessarily reads the same cell voltage.



## 8. Move the input and inject current



Set Click action to Move injection site and click a branch. The white ring moves to that location. Enter a current amplitude and pulse length, then press Inject now.

The pulse begins at the next simulated step. If paused, use Run or Step to advance through it. Pulse length is rounded to the fixed half-millisecond timestep.

Positive and negative pulses are allowed. Larger magnitude usually gives a stronger response, while longer pulses supply current for longer. Repeated or overlapping pulses can interact.

The Inject into selected cell shortcut also allows bridge or balancer stimulation. For an A N N cell it selects the first exported compartment, not an identified soma, so a branch click gives more precise control.



## 9. Read the voltage colors and graph



Voltage means deviation from a resting reference. Zero millivolts means at that reference, not a biological absolute membrane voltage of zero.

In voltage color mode, blue is negative, muted teal is near rest, and gold is positive. The chosen plus-or-minus range controls the display; values beyond it clip to the endpoint color.

The time plot uses milliseconds horizontally and millivolts vertically. All selected probes share one automatically scaled vertical axis, so large responses may visually hide small ones.

The table gives the latest numerical readings. The window selector changes how much recent time is visible. These voltages are not spike rates, calcium fluorescence, ciliary beat frequencies, or swimming predictions.



## 10. Radius changes two things at once



Cable radius defaults to zero point two five micrometers everywhere along the A N N arbors. This uniform radius is assumed, not a measured diameter profile.

Increasing radius increases membrane area, capacitance, and total leak at fixed time constant. It also increases internal cable conductance with the square of radius.

Doubling radius therefore doubles capacitance but quadruples internal conductance. Current can spread more effectively, but the same input is distributed across more membrane.

The resulting voltage is not guaranteed to increase. Radius changes geometry-dependent electrical load and coupling together, so compare several probes and use resets for controlled trials.



## 11. Internal resistance and coupling



Axial resistivity defaults to one hundred fifty ohm centimeters. It represents how difficult it is for current to flow through the internal material.

Increasing resistivity weakens coupling between neighboring compartments. Decreasing it strengthens that coupling and tends to reduce spatial voltage differences.

The internal conductance multiplier defaults to one and changes coupling directly without adding membrane area. Zero removes all internal cable coupling, an extreme mathematical scenario.

With sufficiently strong internal coupling, the spatial model approaches the point model. This limiting behavior is one of the implementation checks, not evidence that the chosen default is physiological.



## 12. Capacitance, leak and persistence



Capacitance describes charge storage per voltage change. Specific capacitance defaults to one microfarad per square centimeter for the A N N membrane.

Bridge and balancer cells each have an assumed capacitance of zero point five picofarads. Increasing capacitance generally reduces voltage response to a fixed current.

The leak time constant defaults to twenty milliseconds. A larger value means weaker leak at fixed capacitance, allowing responses to persist longer.

In this code leak conductance equals capacitance divided by the time constant. Changing capacitance while keeping the time constant fixed also changes leak. You are not varying capacitance alone.



## 13. What a synaptic contact contributes



Each anatomical contact contributes a current determined by its presynaptic voltage. The current scale defaults to zero point zero five picoamperes per contact.

This is a saturating scale, not a constant current delivered at every moment. The model uses a hyperbolic tangent of presynaptic voltage divided by a five millivolt scale.

Increasing current scale strengthens contact effects. Setting it to zero switches off synaptic currents. Increasing the voltage scale makes contacts less sensitive to small voltage changes.

Every contact is assigned equal strength in this release. Multiple contacts add contributions, but anatomical contact count is not a measured functional weight.



## 14. Signs and feedback are hypotheses



The sign setting applies uniformly to every contact. Positive sign makes positive presynaptic deviations generate positive current; negative sign reverses that relationship.

This signed approximation does not identify neurotransmitters or reproduce receptor kinetics. Negative input and feedback loops can make the network response more complicated than simply excitation versus inhibition.

Turning off bridge-to-A N N feedback disables forty four contacts. It does not remove bridge cells or their other connections.

Press Apply parameters now to activate changed settings. Existing voltages are retained and the change is logged. To compare from the same resting state, reset before each trial.



## 15. Cut a cable at a location you choose



Set Click action to Select cable cut and click a branch. A gold cross marks the disconnected internal edge. Selecting a different edge restores the previous one.

You can instead choose the balanced default cut, then press Apply. That preset approximately divides the central arbor into two similarly sized portions of traced cable length.

The cut retains membrane and synaptic contacts on both sides. It models neither cell death nor wound healing, and other cells can still provide indirect paths between the separated regions.

The point comparison cannot represent a local cable cut, so it remains internally equipotential. A local voltage can increase after a cut because there are fewer internal paths for injected current to spread.



## 16. A clean experiment to try



For a first experiment, choose two distant sites on the central A N N and add local probes. Also add that cell mean and one balancer-cell readout.

Reset, inject a pulse at one site, and compare local responses with the cell mean. Then export your results before changing the probe set.

Change one parameter, such as axial resistivity, reset again, and repeat the same input. Keep the sites, pulse amplitude, pulse duration, and other parameters fixed.

Look for differences in peak amplitude, timing, and decay. A reproducible model difference suggests a useful measurement to investigate; it does not yet prove what the real animal does.



## 17. Save traces and the experiment log



Export traces C S V saves the selected probes for both models at each retained half-millisecond time step. The browser retains up to ten thousand recent samples, which is five seconds of simulated time.

Export experiment J S O N saves initial and current settings, interventions, probe-set history, and the retained traces. It is an experiment log, not a resumable full-state checkpoint.

Export before resetting, clearing the window, or changing probes. The plot window does not change the underlying retention limit.

The project also includes original source files, provenance hashes, the Python reference, and numerical tests. Keep these with your exports so the anatomical and mathematical assumptions remain traceable.



## 18. The original wiring and comparison figures



The earlier E X P L O R E page remains in the package as a fixed reference. Its wiring matrix uses rows as senders and columns as receivers.

An entry counts contacts, not neurons, spikes, or measured synaptic strength. That anatomical matrix does not change when you disable feedback in a simulation.

The saved comparison image shows an x y projection of the three arbors, then central A N N mean voltage and balancer mean voltage over time. The time plots use different vertical scales.

The original eighteen-scenario sensitivity table compares radii, signs and current scales. Its root mean squared difference measures disagreement between models, not error against experimental truth or a confidence interval.



## 19. What real time does—and does not—mean



The live solver uses a direct tree solution for passive membrane and internal currents, with synaptic currents evaluated from the previous time step. The live timestep is fixed at zero point five milliseconds.

There is no longer a fixed one hundred sixty millisecond endpoint. The simulation continues until you pause or reset, and switching away from the tab pauses it.

The JavaScript solver was compared with the Python SciPy reference across full compartment states, perturbations, and parameter changes. Their close agreement checks implementation consistency.

Extreme parameters and long runs still require convergence checks and biological constraints. The next scientific goal is a prediction that survives plausible assumptions and can be compared with an independent experiment.

