# Ink/Stitch Satin Tools

Inkscape extensions for building and adjusting **satin columns** for [Ink/Stitch](https://inkstitch.org). A satin column here is one path made of two rail subpaths plus the cross-tie "rung" subpaths that Ink/Stitch uses to shape the stitches. Sometimes it is called a ladder.

| Menu entry | What it does | Extension ID |
|---|---|---|
| **Satin Rail Tools** | One dialog with a tab for seven of the tools below | `maps_edge.satin_tools` |
| **Edit Ladder Rungs** | Builds the rungs for two rails, then adds or removes rungs | `maps_edge.rung_rail_satin` |
| **Add Rung at Node** | Adds one rung at a node you pick on a rail | `maps_edge.add_rung_at_node` |
| **Outset / Inset Satin Column** | Makes a satin column wider or narrower without redrawing it | `maps_edge.outset_satin` |
| **Cut Satin at Nodes** | Splits a satin column into separate columns at nodes, a line, or a shape, or separates a column already in pieces | `maps_edge.cut_satin_at_nodes` |
| **Join Satin** | Joins satin columns that meet end to end into one column | `maps_edge.join_satin` |
| **Satin Center Line** | Draws a line down the middle of a satin column, in its style, just under it | `maps_edge.satin_center_line` |
| **Cross-Splice Rails** | Rejoins 2 or 4 open strokes at a crossing into rails, and builds rungs | `maps_edge.cross_splice_satin` |
| **Normalize Stitch Direction** | Flips rail direction across several columns so they stitch through in one direction | `maps_edge.normalize_stitch_direction` |
| **Weave Paths** | Cuts gaps in crossing strands so they read as woven over and under | `maps_edge.weave_paths` |
| **Touch Rails at Crossing** | Extends the under column's rail tips to touch the over column | `maps_edge.touch_rails` |
| **Taper Join at Crossing** | Tapers the under column's rails to one point on the over column | `maps_edge.taper_join` |
| **Clip Under to Over** | Clips the under column at the over column so satin fill reaches the tip | `maps_edge.clip_to_over` |
| **Reroute Rail Along Shape** | Sends a rail around a closed shape instead of through it, like Difference does for a filled area | `maps_edge.reroute_rail` |
| Angle Rungs | Sets every rung to one angle across the column. No menu entry of its own: it is a tab in **Satin Rail Tools** | (in `maps_edge.satin_tools`) |

All menu entries are under **Extensions > Map's Edge Tools > Embroidery > Satin**. Every tool is also a button on the [Map's Edge Tool Panel](https://github.com/theRoadLessOrdinary/maps-edge-panel).

The tools write the same `inkstitch:satin_column` attribute Ink/Stitch writes itself, so the result is recognized as a satin column right away. Ink/Stitch is only needed to preview and stitch the result; these tools do not call it.

## Requirements

- Inkscape 1.x (developed and used on Inkscape 1.4.4, Flatpak, on Linux; other platforms have not been tried)
- Python 3 with `numpy` available to Inkscape's Python (every tool imports it)
- Ink/Stitch installed, to stitch the results

## Install

Copy every file in this repository into a folder inside Inkscape's user extensions directory, then restart Inkscape. To find that directory, open **Edit > Preferences > System** and look at **User extensions**. Typical locations:

| Install type | User extensions folder |
|---|---|
| Linux (native) | `~/.config/inkscape/extensions/` |
| Linux (Flatpak) | `~/.var/app/org.inkscape.Inkscape/config/inkscape/extensions/` |
| macOS | `~/Library/Application Support/org.inkscape.Inkscape/config/inkscape/extensions/` |
| Windows | `%APPDATA%\inkscape\extensions\` |

For example, on native Linux:

```bash
git clone https://github.com/theRoadLessOrdinary/inkstitch-satin-tools.git \
    ~/.config/inkscape/extensions/inkstitch-satin-tools
```

Keep all the `.py` and `.inx` files together in one folder. To use these tools from the [Map's Edge Tool Panel](https://github.com/theRoadLessOrdinary/maps-edge-panel), name the folder `rung-rail-satin`. The tools import each other (`rung_rail_satin.py` and `bezier_geom.py` are shared code), so copying a single tool by itself will not work. To uninstall, delete the folder and restart Inkscape.

## Typical workflow

1. Draw the two rails of a satin column as open paths.
2. Select both and run **Edit Ladder Rungs** in **Add** mode. Rungs are built and the path is marked as a satin column.
3. Tune rung density with **Add** (more rungs) or **Decimate** (fewer rungs), and add single rungs with **Add Rung at Node**.
4. Set the rung angle with **Angle Rungs** and the width with **Outset / Inset Satin Column**.
5. For a knot or braid: **Weave Paths** on the strands first, **Cross-Splice Rails** where strokes need to swap, then **Clip Under to Over** (or **Touch Rails** / **Taper Join**) where an under piece meets an over piece.
6. Split a column with **Cut Satin at Nodes**.
7. With several columns that will stitch one after another, select them and run **Normalize Stitch Direction**.
8. Open Ink/Stitch's Params and simulator to check the result.

## Which subpaths are the rails

Ink/Stitch does not require the two rails to come first in the path, and after editing they often do not. The tools differ in how they find them:

- **Outset / Inset** uses Ink/Stitch's own rule: the rails are the subpaths that cross more than two other subpaths, or the two longest when that is unclear.
- **Add Rung at Node, Cut Satin at Nodes, Touch Rails, Taper Join, Clip Under to Over, Reroute Rail Along Shape** take the two longest subpaths. Reroute Rail Along Shape treats a path of one or two subpaths as all rails.
- **Edit Ladder Rungs, Angle Rungs, Normalize Stitch Direction** take the first two subpaths.

None of them reads `inkstitch:satin_column`: any path with two or more subpaths is treated as a satin column.

## Satin Rail Tools

One dialog with a tab per tool: Ladder Rungs, Add Rung at Node, Angle Rungs, Cross-Splice Rails, Cut Satin, Normalize Direction, Outset / Inset. Each tab runs the same code as that tool's own menu entry, on the current selection. The separate menu entries still work.

## Edit Ladder Rungs

**Select** one or more ladders (2 rail subpaths plus rungs) and/or separate open rail paths in pairs. Separate rails pair up in the order you Shift+click them: rails 1 and 2 become one ladder, rails 3 and 4 the next, and so on. Each pair is combined into one path: the result keeps the path that is later in the document and the other is deleted. A selected path that already has two rails is edited on its own and does not count toward the pairs. A rubber-band selection or Ctrl+A has no click order, so use Shift+click when selecting separate rails.

Nothing is changed if the number of separate rails is odd, or if a pair's rails are far apart for their length (average gap more than half the longer rail's length), since that usually means a rail was clicked out of order. The message names the rails involved.

**Modes**

- **Add**
  - With only the two rails and no rungs yet, this builds the full set of rungs (see the rules below).
  - With rungs already present, it inserts one new rung midway between every adjacent pair of rungs, which halves the spacing. Run it again to halve it again.
- **Decimate** removes every other rung, keeping the 1st, 3rd, 5th and so on. Run it again to thin further.
- **Fit to rails** trims or extends every rung so it ends 0.125 mm past each rail, keeping the rung's position and angle. The overshoot is the same as new rungs get: 0.125 mm past each rail, on each side. Use it after moving or resizing rails. All rungs are fitted, including the 4 at each end. A rung whose line no longer crosses both rails is left as it is, and the tool says how many.
- **Curls:** where a column curls back on itself, a rung end runs at most halfway across the gap to the next turn, so each rung crosses only its own two rails (Ink/Stitch misreads the width otherwise). This applies to new rungs and to Fit to rails. Leave at least 0.6 mm between turns, about 1.5 times an average Floriani thread, so the ends have room. The tool does not check the gap: getting it wide enough is up to the designer. A rung must also lie inside the column between its two rails and cross no other rail on the way; near a curl's center the tool skips any spot where it would reach across the gap instead, such as two rail tips facing each other across it. At a curl-back, where a rail bends back across a rung's way so the rung meets that same rail a second time before it reaches the other rail, the rung ends at that second crossing: it crosses the one rail twice, which is right there.
- **Remove all** deletes every rung and keeps the two rails, including the rungs at each end that Add and Decimate never touch.
- **Angle** sets every rung to the Degrees setting across the column (formerly the separate Angle Rungs tool, which is now a Ladder Rungs mode). 0 is square across, positive leans clockwise. Rungs inside a taper, and any that would run off the column or cross another rung, are left as they are.
- **Number subpaths** puts a tiny number beside each subpath of the selected paths, counted from 1 in the order they appear in the path's d, in a "Subpath numbers" layer carrying Ink/Stitch's ignore-layer command. The path also gets one attribute per subpath, mapsedge:sb0001 and so on, holding that subpath's coordinates; each number carries mapsedge:of (the path id) and mapsedge:sb. The numbers are a snapshot: number again after editing. Numbering again replaces the old numbers and says if they had gone out of date.
- **Remove numbers** deletes the numbers and attributes from the selected paths, or from every path when nothing is selected, and removes the layer once it is empty.

**Protected ends.** The first 4 and last 4 rungs (in order along the rails) are never touched by Add or Decimate, so tapers and end treatments stay as they were generated. Add needs at least 8 rungs to have any gap to fill; Decimate needs at least 10 rungs before it removes anything. When there is nothing to do, the tool says so and leaves the path unchanged.

Rungs inserted by Add go halfway between their two neighbors on both rails, at the same both-rails pairing the full build uses (see below), and end 0.125 mm past each rail.

**How the first build works.** In the rules below, A is the measured rail-to-rail distance (the median across the middle 70% of the rail) and B is the path's stroke width (0.254 if none is set).

- Body rungs join matching points on the two rails, chosen so each rung is as close to square to **both** rails as it can be. The matching only ever moves forward along each rail, so rungs never cross and never jump to another turn of a rail that curls. On parallel rails this is the same as square across. Where both rails curve around a common center, as in a curl or spiral, the rungs fan out like spokes. Inside a taper, rungs run parallel to the last non-taper rung.
- A **taper** is detected from the path structure: at an end where one rail has more than one segment and its last segment lands at that end, that last segment is the taper. It does not matter whether the segment is straight or curved.
- A taper zone gets 3 rungs, placed no closer than 2B to the convergence point on either rail. A 4th rung is added in the middle only when the taper is clearly longer than A (over 1.15 x A).
- Every rung, taper and end rungs included, ends 0.125 mm past each rail on both sides. Where a column curls back on itself, an end can be shorter (see Curls above).
- An end where the rails meet at a point without a corner (a smooth tip, such as the tail of a curl) is treated like a taper: the body rungs continue right down into the tip, and the last one sits where the column is still at least 4B wide, at least 2B from the point.
- An end with no taper whose rails do not meet gets one final rung parallel to the straight line between the two raw endpoints, pulled back at least 4B from both rails.
- Body rungs are spaced 3A apart on gentle curves and straights, and 1A apart on tighter curves, measured along whichever rail is moving faster at that point (the outside of a curve).
- The two rails are the two longest subpaths, wherever they sit in the path; every other subpath is a rung. Ladder Rungs, Angle Rungs and Normalize Direction put the rails first when they rewrite a path.
- Rails drawn in opposite directions are detected and normalized before rungs are built.

Treat a first build as a draft. The taper classification and the 3-versus-4 rung threshold are heuristics worked out against a small set of real designs, so check the ends.

## Add Rung at Node

In the Node tool (N):

- click one node on either rail where you want a new rung, or
- click two nodes on the same rail to put a rung halfway between them,

then run **Add Rung at Node**. A single-node rung is placed 2 x the stroke width along the rail from the node, square to that rail, crossing to the other rail with an overshoot of 1 x the stroke width.

It has no dialog, so it can be bound to a key in **Edit > Preferences > Interface > Keyboard Shortcuts**.

## Angle Rungs

**Select** one or more satin columns. Every rung is set to the same angle across the column.

- 0 is square across. Positive leans the rungs clockwise on screen, negative counterclockwise (-80 to 80 degrees, default 30).
- The angle is absolute, not added to what is there, so running it again with another value re-aims the rungs.
- Each rung turns about the middle of the column where it sits, and keeps its overshoot past each rail.
- Left as they are: rungs inside a taper (where the rails are more than 20 degrees apart), and rungs whose new line would miss a rail or cross another rung. A message says how many.

## Outset / Inset Satin Column

**Select** one or more satin columns. Every selected path is processed.

**Parameter:** *Distance per rail* in px, from -1000 to 1000 (default 1.0). Positive makes the column wider, negative makes it thinner. The distance applies to each rail, so the column changes width by twice the value. It is measured on the canvas, after any scaling on the path and its groups.

**How it works**

- Each rail moves away from the other rail by the distance, along its own normal. One outward side is used for the whole rail.
- Node handles are adjusted for the curve of each segment, so bends stay even. Where a segment would still be off by more than 0.05 mm, it is split and done again, so rails can gain nodes.
- At a corner the node moves to where both offset sides meet (limited to twice the distance).
- On an inward bend tighter than the distance, the moved rail would fold back on itself; those small loops are cut out.
- Straight rungs are moved onto the new rails and keep the overshoot they had.

**Pointed tips stay pointed.** A rail end within 0.3 mm of the other rail counts as a pointed tip.

- Outset grows from zero at a pointed tip to the full distance over 3 mm.
- Inset never takes more than 40% of the column's width at any point, so the rails never cross. On a column narrower than twice the inset, that stretch gets less than the full inset.
- If the rails would still cross just short of a tip, both are cut where they cross, so they meet there.
- A rung that ends up not reaching both rails (usually right at a tip) is counted in a message, so you can check it or delete it.

Each run is a single Inkscape undo step.

## Cut Satin at Nodes

Splits satin columns into separate satin columns. It works on every satin column in the selection: each column a line or shape crosses is cut, and columns it misses are left alone. Each piece gets its own rail pieces and the rungs on its side of the cut, and keeps the original's style and Ink/Stitch settings. Ways to say where to cut:

- **A line:** select the column(s) plus one or more lines drawn across them. Each line must cross each rail once. The lines are left in place.
- **A closed shape:** draw a rectangle, circle or closed path over the column instead. The column is cut where the outline crosses the rails, the piece inside the shape is deleted with its rungs, and the shape is deleted.
- **Nodes:** in the Node tool, select one node on each rail per cut (nodes can be on several columns at once), or a node on just one rail (the other rail is cut square across from it). Both rails need the same number of nodes, or one rail none. Rail end nodes cannot be cut points.
- **Already in pieces:** select only the column, with no line, shape or nodes. A column whose rails were broken into separate pieces, for example by deleting segments in the Node tool, becomes one satin column per piece, each with its own two rails and the rungs between them. Subpaths that cross belong to the same piece; a piece with no rungs is matched to the nearest other lone rail. Rail pieces that only touch end to end count as separate.

Several cuts in one run are fine. If any column would fail (a line crossing only one of its rails, for example), the message names that column and nothing is changed. A rung that straddles a cut goes to the side it is further into, and a message says how many did.

## Satin Center Line

Select one or more satin columns and run it. Each gets a new path down the
middle of its two rails, with the satin's own style (stroke color, width
and the rest), placed just under the satin in the stacking order. The satin
column is kept. The new path is a plain path, not a satin column, so
Ink/Stitch stitches it as a running stitch. Rungs decide which point on one
rail goes with which on the other. The line has a node wherever a rail or a
rung has one, plus extra nodes only where a curve needs them to stay within
0.05 mm of the true middle.

## Join Satin

The reverse of Cut Satin. Select two or more satin columns that meet end to end, so that each rail's end node sits on, or within 0.5 mm of, a rail end of the next column. Run **Join Satin** and they become one satin column:

- each rail joins the rail on its own side of the next column, and the two end nodes at each joint become one node at their midpoint, keeping both handles;
- every rung is kept;
- the columns can be selected in any order and drawn in any direction. Rails that run opposite ways, or swap sides from one column to the next, are lined up first;
- the joined column keeps the style and Ink/Stitch settings of the lowest selected column and takes its place. The others are removed.

Columns in a closed loop become one column that starts and ends at the same joint. If a column doesn't meet the others, a message names it and nothing changes.
With nodes selected in the Node tool, only the joints where a selected node
sits are joined; ends that touch elsewhere stay apart. Two pieces of a loop
joined at one joint make one column that starts and ends at the other.


## Cross-Splice Rails

For a crossing (a Celtic knot, a braid) where the strokes need to swap which one continues into which.

**Select** exactly 2 or exactly 4 open strokes around **one** crossing. Each must be a single open subpath. Shapes that can be converted to paths are converted.

- **4 strokes:** the tool tries the 3 ways to split them into two pairs and uses the split with the smallest total end-to-end distance. Each pair is spliced into one continuous rail, giving 2 rails, and rungs are generated in the same step (as **Edit Ladder Rungs** would).
- **2 strokes:** they are spliced into one open path and no rungs are built. Use this to feed a hand-assembled rail into **Edit Ladder Rungs**.

**Pairing is by distance, not by direction.** An earlier version scored pairs by tangent alignment and matched strokes across a whole knot design to the wrong loops, because the design repeats by rotation. Nearest-distance matching does not have that problem when the selection is limited to one crossing, which is why you should select only the strokes at that crossing.

**The joint.** A genuine taper at the joining end is dropped entirely (not bridged over) and replaced by one straight line to the new joint. A non-tapered end is extended by a straight segment. The old curve's arc is not matched at the joint. If you want it smoother, add nodes a couple of millimeters either side of the joint and delete the joint nodes so Inkscape's curve fitting smooths it. Strokes that carry different transforms are handled: each stroke's full transform is applied to its path data before splicing.

## Normalize Stitch Direction

**Select** two or more already-built ladder objects. They are processed in document (z-) order, so arrange them in the order they will be stitched.

- The first object is oriented bottom to top: its rails start at the end with the larger Y value (SVG Y grows downward).
- Each following object has both rails reversed together, if needed, so its start is as close as possible to where the previous object's rail ended.

The result is one continuous direction through the whole selection, without jumping to a far point and doubling back. Only the two rail subpaths are ever reversed. Rungs are left alone, since their geometry does not depend on direction. Each selected object must already have its rails and rungs (2 or more subpaths).

## Weave Paths

**Select** two or more open paths that cross each other, such as the strands of a knot or braid. Each strand is cut with a gap wherever it passes under another, alternating over and under along its own length, so the strands read as woven together.

- Strands are handled in selection order. Along each strand, its crossings alternate over, under, over, starting with over. A crossing already decided by an earlier strand is kept.
- A gap is 1.5 x the stroke width of the strand on top.
- Each cut strand becomes separate open paths, one per piece.

## Touch Rails at Crossing

**Select** two satin columns meeting at one crossing: the under piece, which has a cut end, and the over piece it should touch.

- The selected object lower in the stacking order is the under piece, and is the only one changed.
- The under piece's end nearest the over piece is the one fixed. Each of its two rail tips is extended along its own direction (or pulled back) until it touches the over piece's near rail. The end is chosen per rail, so rails drawn in opposite directions work.
- Rungs are left as they are.

Run it once per crossing.

## Taper Join at Crossing

**Select** two satin columns meeting at one crossing, as for Touch Rails. The under piece is the one with a rail tip nearest the other piece.

- If the columns cross at more than about 5 degrees, the under piece's two rails taper together to one shared point on the over piece, with one or two new rungs near the tip. The rungs run parallel to the nearest rail of the first selected path.
- At a shallower graze it does what Touch Rails does.

## Clip Under to Over

**Select** two satin columns meeting at one crossing. The under piece is the one with a rail tip nearest the other piece.

- Each of the under piece's rails is clipped exactly where it enters the over piece.
- At an angled crossing the two rails enter at different points, which would leave the tip without a second rail for Ink/Stitch to fill against. So the shorter rail is extended along the over piece's near rail until it lines up with the longer rail's end, and satin fill reaches all the way to the tip.
- Under-piece rungs with a node inside the over piece are removed.

## Reroute Rail Along Shape

**Select** one closed shape (a path, rectangle, circle or ellipse) and one or more rails: open paths, satin columns, or both. No dialog opens.

- Where a rail runs through the shape, the part inside is replaced by the shape's edge between the two points where the rail crosses it. The rail now goes around the shape instead of through it. This is what Path > Difference does to a filled area, done on an open rail.
- The edge is followed the shorter way around the shape.
- Crossing points are exact: curves are split, not approximated.
- A rail that starts or ends inside the shape is trimmed at the edge.
- A rail that never crosses the shape's edge is left as it is.
- On a satin column only the two rails are rerouted. Rungs are left as they are, so a rung inside the shape stays there.
- The shape is deleted, as Difference does.

## Files

| File | Purpose |
|---|---|
| `satin_tools.inx`, `satin_tools.py` | **Satin Rail Tools** dialog; runs the other tools' code |
| `rung_rail_satin.inx`, `edit_ladder_rungs.py` | **Edit Ladder Rungs** |
| `rung_rail_satin.py` | Shared module: rung generation, rail geometry, taper classification. Not a menu entry by itself |
| `bezier_geom.py` | Shared Bezier helpers: intersection, flattening, splitting. Uses only `numpy` |
| `add_rung_at_node.inx`, `add_rung_at_node.py` | **Add Rung at Node** |
| `angle_rungs.py` | **Angle Rungs** (reached through Satin Rail Tools) |
| `outset_satin.inx`, `outset_satin.py` | **Outset / Inset Satin Column** |
| `cut_satin_at_nodes.inx`, `cut_satin_at_nodes.py` | **Cut Satin at Nodes** |
| `join_satin.inx`, `join_satin.py` | **Join Satin** |
| `satin_center_line.inx`, `satin_center_line.py` | **Satin Center Line** |
| `cross_splice_satin.inx`, `cross_splice_satin.py` | **Cross-Splice Rails** |
| `normalize_stitch_direction.inx`, `normalize_stitch_direction.py` | **Normalize Stitch Direction** |
| `weave_paths.inx`, `weave_paths.py` | **Weave Paths** |
| `touch_rails.inx`, `touch_rails.py` | **Touch Rails at Crossing** |
| `taper_join.inx`, `taper_join.py` | **Taper Join at Crossing** |
| `clip_to_over.inx`, `clip_to_over.py` | **Clip Under to Over** |
| `reroute_rail.inx`, `reroute_rail.py` | **Reroute Rail Along Shape** |

## Known limits

- Rung generation is a rule set worked out by hand against real designs. It has not been tested beyond those, so review the ends of every column.
- Outset / Inset adds nodes to rails where needed for accuracy (about twice as many on one test design).
- Cross-Splice joins with straight connectors and does not match curve arcs at the joint.
- Touch Rails extends rail tips with straight lines.
- Reroute Rail Along Shape always takes the shorter way around the shape. When the rail enters and leaves through different holes of a shape, that piece is cut out instead of rerouted.
- Weave Paths and Taper Join rely on selection order, which Inkscape does not guarantee in every case.
- Only the stated selections are supported. Anything else stops with an error message and leaves the document unchanged.

## License

GPL-3.0. See [LICENSE](LICENSE).
