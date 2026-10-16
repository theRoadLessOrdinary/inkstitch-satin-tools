# Ink/Stitch Satin Tools

Four Inkscape extensions for building and adjusting **satin column ladders** for [Ink/Stitch](https://inkstitch.org). A ladder is one path made of two rail subpaths plus the cross-tie "rung" subpaths that Ink/Stitch uses to shape a satin column.

| Menu entry | What it does | Extension ID |
|---|---|---|
| **Edit Ladder Rungs** | Builds the rungs for two rails, then adds or removes rungs on an existing ladder | `maps_edge.rung_rail_satin` |
| **Outset / Inset Satin Column** | Makes a built satin column wider or narrower without redrawing it | `maps_edge.outset_satin` |
| **Cross-Splice Rails** | Rejoins 2 or 4 open strokes at a crossing into rails, and builds rungs | `maps_edge.cross_splice_satin` |
| **Normalize Stitch Direction** | Flips rail direction across several ladders so they stitch through in one continuous direction | `maps_edge.normalize_stitch_direction` |

All four appear under **Extensions > Map's Edge Tools > Embroidery > Satin**. They write the same `inkstitch:satin_column` attributes Ink/Stitch writes itself, so the result is recognized as a satin column right away. Ink/Stitch is only needed to preview and stitch the result; these tools do not call it.

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

Keep all the `.py` and `.inx` files together in one folder. The tools import each other (`rung_rail_satin.py` and `bezier_geom.py` are shared code), so copying a single tool by itself will not work. To uninstall, delete the folder and restart Inkscape.

## Typical workflow

1. Draw the two rails of a satin column as open paths.
2. Select both and run **Edit Ladder Rungs** in **Add** mode. Rungs are built and the path is marked as a satin column.
3. Tune rung density with **Add** (more rungs) or **Decimate** (fewer rungs), as many times as needed.
4. Adjust the width with **Outset / Inset Satin Column**.
5. Where strokes cross (a knot or braid), use **Cross-Splice Rails** to reconnect them into rails first.
6. With several ladders that will stitch one after another, select them and run **Normalize Stitch Direction**.
7. Open Ink/Stitch's Params and simulator to check the result.

## Edit Ladder Rungs

**Select** either one ladder (2 rail subpaths plus rungs) or two separate open rail paths. Two separate paths are combined into one: the result keeps the path that is later in the document and the other is deleted.

**Modes**

- **Add**
  - With only the two rails and no rungs yet, this builds the full set of rungs (see the rules below).
  - With rungs already present, it inserts one new rung midway between every adjacent pair of rungs, which halves the spacing. Run it again to halve it again.
- **Decimate** removes every other rung, keeping the 1st, 3rd, 5th and so on. Run it again to thin further.

**Protected ends.** The first 4 and last 4 rungs (in order along the rails) are never touched by Add or Decimate, so tapers and end treatments stay as they were generated. Add needs at least 8 rungs to have any gap to fill; Decimate needs at least 10 rungs before it removes anything. When there is nothing to do, the tool says so and leaves the path unchanged.

Rungs inserted by Add run perpendicular to the first rail at the midpoint and overshoot each rail by 1 x the stroke width.

**How the first build works.** In the rules below, A is the measured rail-to-rail distance (the median across the middle 70% of the rail) and B is the path's stroke width (0.254 if none is set).

- Rungs are perpendicular to the local rail tangent, except inside a taper, where they run parallel to the last non-taper rung.
- A **taper** is detected from the path structure: at an end where one rail has more than one segment and its last segment lands at that end, that last segment is the taper. It does not matter whether the segment is straight or curved.
- A taper zone gets 3 rungs, placed no closer than 2B to the convergence point on either rail. A 4th rung is added in the middle only when the taper is clearly longer than A (over 1.15 x A).
- Taper rungs overshoot by 1B on the non-taper rail and 8B on the taper rail, so the node marker there does not hide them.
- An end with no taper gets one final rung parallel to the straight line between the two raw endpoints, pulled back at least 4B from both rails, overshooting 4B on both sides.
- Interior rungs are spaced 2A apart on gentle curves and straights, and 1A apart on tighter curves. Default overshoot is 1B on both sides.
- Rails drawn in opposite directions are detected and normalized before rungs are built.

Treat a first build as a draft. The taper classification and the 3-versus-4 rung threshold are heuristics worked out against a small set of real designs, so check the ends.

## Outset / Inset Satin Column

**Select** one or more already-built ladder objects (2 rail subpaths followed by rung subpaths). Every selected path is processed.

**Parameter:** *Distance per rail* in px, from -1000 to 1000 (default 1.0). Positive makes the column wider, negative makes it thinner. The distance applies to each rail, so the finished column changes width by twice the value.

**How it works**

- Each rail is pushed away from the other rail along its local normal. The direction is worked out point by point, so it is correct when the rails curve or change spacing.
- Each rung endpoint moves by the same amount as the nearest rail point, so rungs stay attached to the rails and keep their overshoot.
- The distance is corrected for the object's own scale transform, so 1.0 is 1 px in document units even on a scaled object.
- Nodes and their handles move as a group. This is an approximation of a true offset curve, not an exact one. On very tight curves a large distance can distort the shape, so use several small steps and check the result.

Each run is a single Inkscape undo step.

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

## Files

| File | Purpose |
|---|---|
| `rung_rail_satin.inx` | Menu definition for **Edit Ladder Rungs** (runs `edit_ladder_rungs.py`) |
| `edit_ladder_rungs.py` | Add / Decimate logic |
| `rung_rail_satin.py` | Shared module: rung generation, rail geometry, taper classification. Not a menu entry by itself |
| `bezier_geom.py` | Shared Bezier helpers: intersection, flattening, splitting. Uses only `numpy` |
| `outset_satin.inx`, `outset_satin.py` | **Outset / Inset Satin Column** |
| `cross_splice_satin.inx`, `cross_splice_satin.py` | **Cross-Splice Rails** |
| `normalize_stitch_direction.inx`, `normalize_stitch_direction.py` | **Normalize Stitch Direction** |

## Known limits

- Rung generation is a rule set worked out by hand against real designs. It has not been tested beyond those, so review the ends of every ladder.
- The outset/inset offset is per-node, not an exact offset curve.
- Cross-Splice joins with straight connectors and does not match curve arcs at the joint.
- Only the stated selections are supported. Anything else stops with an error message and leaves the document unchanged.

## License

GPL-3.0. See [LICENSE](LICENSE).
