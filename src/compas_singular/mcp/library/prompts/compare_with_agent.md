Run the comparison between the two meshing agents on the mesh in hand.

There are two independent implementations in this repository:

- `compas_singular.agent` -- an in-process tool loop that calls a model itself.
- `compas_singular.mcp` -- this server, driven by whichever client you are.

They share the library underneath and share nothing else. The point of the
comparison is to find out whether the difference in arrangement shows up in the
result.

For the mesh in hand:

1. `rhino_pull` (or load the same mesh from disk), then `inspect` and record the
   starting triple: min angle, max angle, max aspect.
2. Improve it using only this server's tools. Record every step.
3. `inspect` again and record the finishing triple and the step count.
4. `save_example` with an honest verdict and a remark saying what the deciding
   move was.

Then report a table: starting triple, finishing triple, number of steps, and
which tool produced the largest single improvement. Do not push to Rhino -- this
is a measurement, not an edit.

State plainly if the mesh could not be improved. A run that changes nothing is a
valid result and belongs in the corpus with verdict `acceptable` and a remark
saying why.
