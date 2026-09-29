# Data audit

GitHub tables and scripts are pinned to commit
`ccdc5b2d8e897dcb4d1df1d4a651b6925844320a` in the Jékely lab repository linked in
README. Three public CATMAID project 33 arbors were retrieved separately on
2026-09-27. Exact files, hashes and URLs are in data/provenance.json.

The live arbors and pinned tables are not guaranteed to be identical versions.
All modeled ANN synaptic treenode IDs were checked against the downloaded arbors.
The bundled snapshot and hashes define this release; reruns do not query servers.

| Skeleton ID | Older source name | Display name |
|---|---|---|
| 4019221 | SSN_Q1Q2Q3Q4 | ANN Q1-4 |
| 4018492 | SSN_Q1Q2 | ANN Q1Q2 |
| 4018688 | SSN_Q3Q4 | ANN Q3Q4 |

The mapping is established by shared connector IDs between stats_synapse.csv
and Figure3_source_data1.csv and by Figure3.R. The source label SSN here must not
be confused with the separate SNN label.

The raw endpoint table has 1,132 records and 430 connector IDs: 396 presynaptic
records and 736 postsynaptic records. Joining connectors with both sides yields
684 directed contacts over all cell types. Thirty-four connectors have only
postsynaptic endpoints (52 records); one has only a presynaptic endpoint. These
are reported in results/audit.json, never assigned guessed partners.

The selected circuit contains 3 ANN, 26 bridge and 129 balancer cells. It has
280 contacts over 113 directed pairs, and no contacts leaving balancer cells.
The Figure 3 grouped source matrix also sums to 280. Figure3.R establishes that
Var1 sends and Var2 receives. The matrix includes a bridge-to-balancer contact
and a local ANN autapse; these are preserved despite simpler verbal descriptions.
Matching totals does not establish every group assignment; a full quadrant-label
audit is a next step.

The retained older catmaid-connectivity-matrix.json contains different skeleton
IDs and is not used to define the model. Live CATMAID connector entries are also
not mixed into the pinned contact table.

| Arbor | Raw nodes | Traced length in original audit (µm) |
|---|---:|---:|
| ANN Q1-4 | 8,280 | 2,324.082 |
| ANN Q1Q2 | 5,042 | 1,349.856 |
| ANN Q3Q4 | 3,127 | 847.417 |

After merging three zero-length edges, there are 16,446 cable compartments and
155 other cell compartments. Each export is a rooted tree and includes uncertain
end/continuation tags. The tree does not establish all anastomotic loops. Most
source radii are unspecified (-1); a uniform model radius is declared explicitly.
No smoothing, proximity-based joining or invented missing synapses are used.

Figure 4 source files contain ciliary measurements, not neural voltage. Before
fitting, resolve larva IDs, repeated recordings, anatomical planes, elapsed-time
calibration and preprocessing. Raw frame indices alone do not establish seconds.
