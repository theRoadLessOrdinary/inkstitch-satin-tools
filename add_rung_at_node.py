#!/usr/bin/env python3
"""
Add Rung at Node (Ink/Stitch satin column helper)

Switch to the Node tool (N), then either:
 - click a single node on ONE of the ladder's two rails at the position
   where you want a new rung, or
 - click two nodes on the SAME rail to place the rung equidistant between
   them (midpoint arc-length position, no 2*B offset applied),
then run this extension (Extensions menu, or a keyboard shortcut bound to
it in Preferences -> Interface -> Keyboard Shortcuts - about as fast as a
right-click would be, which Inkscape's extension system has no hook for:
INX menu placement only supports the Extensions menu, confirmed against
inkex's own docs).

Inkscape passes the id/subpath/node of whatever node is selected in the Node
tool to every extension automatically via --selected-nodes=id:subpath:node
(0-based), regardless of what the .inx declares - confirmed straight from
Inkscape's own C++ source (src/extension/implementation/script.cpp). Since
the clicked node's index tells us exactly which rail segment boundary it
sits on, its arc-length position along the rail is exact (Rail.seg_bounds),
not a nearest-point approximation.

New rung = perpendicular to the clicked rail's local tangent, placed 2*B
(stroke width) from the clicked node along the rail (toward whichever end
is nearer, flipping direction if that would run off the rail), crossing to
the other rail, with the same 1*B overshoot convention used for interior
rungs in rung_rail_satin.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, PathElement

from rung_rail_satin import Rail, make_csp, unit_normal, stroke_width, STROKE_DEFAULT_MULT


def local_tangent(rail, s, eps=None):
    """Tangent direction at arc-length s, via a symmetric secant rather than
    the analytic derivative - point_tangent()'s derivative is exactly zero
    at a node joining two straight (degenerate-cubic) segments, since both
    control points collapse onto their endpoints there (same issue flagged
    in the older embroidery-connectors code). A small secant around s sidesteps
    it and still tracks real curvature closely enough for this purpose."""
    if eps is None:
        eps = max(1e-6, 0.02 * rail.length)
    lo = max(0.0, s - eps)
    hi = min(rail.length, s + eps)
    if hi - lo < 1e-9:
        return np.array([1.0, 0.0])
    p_lo, _ = rail.point_tangent(lo)
    p_hi, _ = rail.point_tangent(hi)
    d = np.array(p_hi) - np.array(p_lo)
    n = np.linalg.norm(d)
    return d / n if n > 1e-9 else np.array([1.0, 0.0])


def _s_of_nearest(rail, q):
    """Arc-length position on `rail` nearest to point q. Used as a fallback
    anchor when no perpendicular crossing exists (same approach as
    edit_ladder_rungs.py's s_of_nearest)."""
    best_s, best_d = 0.0, float("inf")
    for i, (ts, arclen, pts) in enumerate(rail._tables):
        d = np.linalg.norm(pts - q, axis=1)
        j = int(np.argmin(d))
        if d[j] < best_d:
            best_d, best_s = d[j], rail.seg_bounds[i] + arclen[j]
    return best_s


# Shared with the other satin tools; kept importable from here (Cut Satin uses it).
from rung_rail_satin import find_rail_indices


def _validate_rail(csp, sp):
    """Common validation shared by the single-node and two-node entry
    points. Returns other_idx (the companion rail's subpath index)."""
    if len(csp) < 2:
        raise ValueError("This path has fewer than 2 rail subpaths - it isn't a built ladder yet.")

    rail_indices = find_rail_indices(csp)
    if sp not in rail_indices:
        raise ValueError(
            "The clicked node is on subpath "
            f"{sp}, not one of the two rails (identified by arc length as "
            f"subpaths {rail_indices[0]} and {rail_indices[1]}). Rungs "
            "themselves aren't rails - click a node on a rail instead."
        )
    return rail_indices[0] if sp == rail_indices[1] else rail_indices[1]


def _rung_from_position(csp, sp, other_idx, s, _fallback=False):
    """Build a perpendicular rung from arc-length position s on subpath sp,
    crossing to subpath other_idx. Returns (pA, pB, u) or raises ValueError."""
    active_rail = Rail(csp[sp])
    other_rail = Rail(csp[other_idx])

    pA, _ = active_rail.point_tangent(s)
    tanA = local_tangent(active_rail, s)
    d = unit_normal(tanA)

    res = other_rail.line_cross(np.array(pA), d)
    if res is None:
        if _fallback:
            raise ValueError(
                "Couldn't find a crossing point from either rail near this "
                "position. Try different nodes."
            )
        # No perpendicular crossing (e.g. the active rail runs past where
        # its shorter companion ends) - snap to the nearest node on the
        # other rail and redo placement anchored there instead, using the
        # exact same perpendicular-rung logic, just with the active rail
        # switched.
        s_b = _s_of_nearest(other_rail, np.array(pA))
        return _rung_from_position(csp, other_idx, sp, s_b, _fallback=True)
    pB, _ = res

    dvec = np.array(pB) - np.array(pA)
    length = np.linalg.norm(dvec)
    if length < 1e-9:
        raise ValueError("The two rails meet exactly at this point - no room for a rung.")
    u = dvec / length

    return pA, pB, u


def compute_rung(csp, sp, nl, B):
    """csp: list of subpaths (as from el.path.to_absolute().to_superpath()).
    sp, nl: 0-based subpath/node index of the clicked node.
    B: stroke width - the rung is placed 2*B away from the clicked node
    along the rail, not exactly at it.
    Returns (endA, endB) world-less/local-space points for the new rung, or
    raises ValueError with a user-facing message."""
    other_idx = _validate_rail(csp, sp)
    active_rail = Rail(csp[sp])

    if nl >= len(active_rail.seg_bounds):
        raise ValueError("That node index doesn't fall on this rail's segment boundaries.")

    s = active_rail.seg_bounds[nl]

    # Place the rung 2*B from the clicked node, toward whichever rail end
    # is nearer; flip direction if that would run off the rail.
    offset = 2.0 * B
    dist_to_start = s
    dist_to_end = active_rail.length - s
    if dist_to_start <= dist_to_end:
        candidate = s - offset
        if candidate < 0.0:
            candidate = s + offset
    else:
        candidate = s + offset
        if candidate > active_rail.length:
            candidate = s - offset
    s = min(max(candidate, 0.0), active_rail.length)

    return _rung_from_position(csp, sp, other_idx, s)


def compute_rung_midpoint(csp, sp, nl_a, nl_b):
    """Two nodes clicked on the same rail: place the rung equidistant
    between them (the midpoint arc-length position), instead of offset
    from a single clicked node."""
    other_idx = _validate_rail(csp, sp)
    active_rail = Rail(csp[sp])

    for nl in (nl_a, nl_b):
        if nl >= len(active_rail.seg_bounds):
            raise ValueError("That node index doesn't fall on this rail's segment boundaries.")

    s = (active_rail.seg_bounds[nl_a] + active_rail.seg_bounds[nl_b]) / 2.0

    return _rung_from_position(csp, sp, other_idx, s)


class AddRungAtNode(EffectExtension):
    def effect(self):
        selected = self.options.selected_nodes
        if len(selected) not in (1, 2):
            raise inkex.AbortExtension(
                "Switch to the Node tool (N), then either click one node "
                "on one of the ladder's two rails at the position you "
                "want the new rung, or click two nodes on the same rail "
                "to place the rung equidistant between them, then run "
                "this again."
            )

        try:
            parsed = [sel.split(":", 2) for sel in selected]
            ids = [p[0] for p in parsed]
            sps = [int(p[1]) for p in parsed]
            nls = [int(p[2]) for p in parsed]
        except (ValueError, IndexError):
            raise inkex.AbortExtension(f"Couldn't parse selected node(s) {selected!r}.")

        if len(selected) == 2 and (ids[0] != ids[1] or sps[0] != sps[1]):
            raise inkex.AbortExtension(
                "Both selected nodes must be on the same rail to place a "
                "rung equidistant between them."
            )

        el = self.svg.getElementById(ids[0])
        if el is None or not isinstance(el, PathElement):
            raise inkex.AbortExtension("The selected node's object isn't a path.")

        csp = [list(subpath) for subpath in el.path.to_absolute().to_superpath()]
        B = stroke_width(el)

        try:
            if len(selected) == 2:
                pA, pB, u = compute_rung_midpoint(csp, sps[0], nls[0], nls[1])
            else:
                pA, pB, u = compute_rung(csp, sps[0], nls[0], B)
        except ValueError as e:
            raise inkex.AbortExtension(str(e))

        end_a = np.array(pA) - u * (STROKE_DEFAULT_MULT * B)
        end_b = np.array(pB) + u * (STROKE_DEFAULT_MULT * B)
        new_rung = [[list(end_a)] * 3, [list(end_b)] * 3]

        el.path = make_csp(csp + [new_rung]).to_path()


if __name__ == "__main__":
    AddRungAtNode().run()
