#!/usr/bin/env python3
"""
Cut Satin at Nodes (Ink/Stitch satin column helper)

Either select a satin column (2 rails + rungs, one combined path) plus
one or more lines drawn across it where you want to cut, or, in the Node
tool (N), select the cut node(s) on the column's rails. Then run this. Each rail is cut at
its selected node(s), and the column is split into separate satin
columns, each made of its own rail pieces plus the rungs that fall on
its side of the cut. Same as doing Break Apart, then re-combining each
side by hand.

Cut line(s): each rail is cut wherever a line crosses it, so each line
has to cross each rail exactly once. Any single-subpath path in the
selection counts as a cut line. The lines are left in place afterwards.

Closed shape(s): draw a closed shape (rectangle, circle, ellipse, polygon
or closed path) over the column instead of a line. The column is cut
wherever the shape's outline crosses the rails, the column piece(s)
inside the shape are deleted along with their rungs, and the shape
itself is deleted. Lines and shapes can be mixed in one run.

Node selection:
 - one node on EACH rail per cut: both rails are cut exactly there, or
 - one node on ONE rail per cut: the other rail is cut where a line
   perpendicular to the first rail at that node crosses it.
Several cuts at once are fine (e.g. 2 nodes per rail = 3 pieces), as long
as both rails get the same number of nodes, or one rail gets none.
Rail end nodes can't be cut points.

Each rung is assigned to a piece by where it crosses each rail. A rung
that crosses rail A on one side of a cut and rail B on the other side
(it straddles the cut line) goes to whichever side it's further into,
and a warning names how many did that.

The pieces keep the original's style and all its attributes (Ink/Stitch
params included); the original path is removed.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, PathElement, Rectangle, Circle, Ellipse, Polygon

import bezier_geom as bg
from rung_rail_satin import Rail, make_csp, unit_normal
from add_rung_at_node import find_rail_indices, local_tangent

# Distance (doc units) within which a computed cut point counts as landing
# on an existing node, so no near-zero-length segment gets created.
NODE_SNAP = 1e-6
# Same idea in distance terms: arc-length lookups are sampled, so a cut
# computed to land exactly on a node can come back a hair off it.
NODE_SNAP_DIST = 1e-3


def split_subpath_at_node(subpath, nl):
    """Split a CSP subpath at node index nl. Both halves share that node."""
    return subpath[:nl + 1], subpath[nl:]


def split_subpath_at_s(rail, subpath, s):
    """Split a CSP subpath at arc length s, inserting a new node there
    (de Casteljau, so the curve shape doesn't change)."""
    i, local_s = rail.seg_index_at(s)
    ts, arclen, _ = rail._tables[i]
    t = float(np.interp(local_s, arclen, ts))
    p0, p1, p2, p3 = rail.segs[i]
    pt = np.array(bg.split_point(p0, p1, p2, p3, t))
    if t < NODE_SNAP or np.linalg.norm(pt - p0) < NODE_SNAP_DIST:
        return split_subpath_at_node(subpath, i)
    if t > 1 - NODE_SNAP or np.linalg.norm(pt - p3) < NODE_SNAP_DIST:
        return split_subpath_at_node(subpath, i + 1)
    left, right = bg.split_cubic(p0, p1, p2, p3, t)
    first = [list(map(list, n)) for n in subpath[:i + 1]]
    second = [list(map(list, n)) for n in subpath[i + 1:]]
    first[-1][2] = list(left[1])
    mid = list(left[3])
    first.append([list(left[2]), mid, list(mid)])
    second[0][0] = list(right[2])
    second.insert(0, [list(mid), list(mid), list(right[1])])
    return first, second


def split_rail(rail, subpath, cuts):
    """cuts: list of (s, node_index_or_None), sorted by s. Returns the
    rail's pieces in order from its start."""
    pieces = []
    rest = subpath
    rest_rail = rail
    consumed_s = 0.0
    consumed_nodes = 0
    for s, nl in cuts:
        if nl is not None:
            first, rest = split_subpath_at_node(rest, nl - consumed_nodes)
        else:
            first, rest = split_subpath_at_s(rest_rail, rest, s - consumed_s)
        pieces.append(first)
        consumed_nodes += len(first) - 1
        consumed_s = s
        rest_rail = Rail(rest)
    pieces.append(rest)
    return pieces


def nearest_s_to_polyline(rail, poly):
    """Arc-length position on `rail` closest to a polyline (a rung). Rungs
    overshoot the rails, so this is the rung's crossing point on that rail."""
    best_s, best_d = 0.0, float("inf")
    for i, (ts, arclen, pts) in enumerate(rail._tables):
        for a, b in zip(poly[:-1], poly[1:]):
            ab = b - a
            denom = float(np.dot(ab, ab))
            if denom < 1e-18:
                proj = np.tile(a, (len(pts), 1))
            else:
                u = np.clip((pts - a) @ ab / denom, 0.0, 1.0)
                proj = a + np.outer(u, ab)
            d = np.linalg.norm(pts - proj, axis=1)
            j = int(np.argmin(d))
            if d[j] < best_d:
                best_d, best_s = d[j], rail.seg_bounds[i] + arclen[j]
    return best_s


def piece_index(s, cut_s):
    return sum(1 for c in cut_s if s > c)


def distance_into_piece(s, cut_s, idx, length):
    lo = cut_s[idx - 1] if idx > 0 else 0.0
    hi = cut_s[idx] if idx < len(cut_s) else length
    return min(s - lo, hi - s)


def rail_setup(el):
    csp = [[list(map(list, n)) for n in sp] for sp in el.path.to_absolute().to_superpath()]
    if len(csp) < 2:
        raise inkex.AbortExtension("This path has fewer than 2 subpaths, so it isn't a satin column.")
    ia, ib = find_rail_indices(csp)
    return csp, ia, ib, {ia: Rail(csp[ia]), ib: Rail(csp[ib])}


def line_crossings(rail, line_csp):
    """Arc-length positions where any segment of the cut line crosses
    `rail`, sorted, with near-duplicates (a hit right on a node, found
    from both neighboring segments) merged."""
    found = []
    for sp in line_csp:
        for j in range(len(sp) - 1):
            c1 = (tuple(sp[j][1]), tuple(sp[j][2]), tuple(sp[j + 1][0]), tuple(sp[j + 1][1]))
            for i, seg in enumerate(rail.segs):
                c2 = tuple(tuple(p) for p in seg)
                for _t1, t2, _pt in bg.dedup_intersections(bg.cubic_intersect(c1, c2)):
                    ts, arclen, _ = rail._tables[i]
                    found.append(rail.seg_bounds[i] + float(np.interp(t2, ts, arclen)))
    found.sort()
    merged = []
    for s in found:
        if not merged or s - merged[-1] > NODE_SNAP_DIST:
            merged.append(s)
    return merged


def is_closed_cutter(el):
    """A closed shape deletes what's inside it; an open line only cuts."""
    if isinstance(el, (Rectangle, Circle, Ellipse, Polygon)):
        return True
    path = el.path.to_absolute()
    if any(seg.letter in "zZ" for seg in path):
        return True
    sp = path.to_superpath()[0]
    return len(sp) > 2 and np.linalg.norm(np.array(sp[0][1]) - np.array(sp[-1][1])) < NODE_SNAP


def close_subpath(sp):
    """Make sure a closed shape's subpath includes its closing segment."""
    sp = [list(map(list, n)) for n in sp]
    if np.linalg.norm(np.array(sp[0][1]) - np.array(sp[-1][1])) >= NODE_SNAP:
        first = list(sp[0][1])
        sp[-1][2] = list(sp[-1][1])
        sp.append([list(first), list(first), list(first)])
    return sp


def shape_polygon(sp):
    poly = []
    for j in range(len(sp) - 1):
        poly.extend(bg.flatten_cubic(sp[j][1], sp[j][2], sp[j + 1][0], sp[j + 1][1])[:-1])
    return poly


class CutSatinAtNodes(EffectExtension):
    def effect(self):
        lines = [e for e in self.svg.selection.filter(PathElement)
                 if len(e.path.to_superpath()) == 1]
        lines += list(self.svg.selection.filter(Rectangle, Circle, Ellipse, Polygon))
        shapes = [e for e in lines if is_closed_cutter(e)]
        polygons = []
        if lines:
            el, csp, ia, ib, rails, cuts, polygons = self.cuts_from_lines(lines, shapes)
        else:
            el, csp, ia, ib, rails, cuts = self.cuts_from_nodes()

        # Rails are often drawn anti-parallel; if so, rail B's order runs
        # backwards relative to rail A's, and so do its pieces.
        sa, ea = np.array(rails[ia].start()), np.array(rails[ia].end())
        sb, eb = np.array(rails[ib].start()), np.array(rails[ib].end())
        reversed_b = (np.linalg.norm(sa - eb) + np.linalg.norm(ea - sb)) < \
                     (np.linalg.norm(sa - sb) + np.linalg.norm(ea - eb))

        pieces_a = split_rail(rails[ia], csp[ia], cuts[ia])
        pieces_b = split_rail(rails[ib], csp[ib], cuts[ib])
        if reversed_b:
            pieces_b.reverse()
        n_pieces = len(pieces_a)

        cut_s_a = [s for s, _ in cuts[ia]]
        cut_s_b = [s for s, _ in cuts[ib]]
        len_b = rails[ib].length

        groups = [[] for _ in range(n_pieces)]
        straddled = 0
        for k, sp in enumerate(csp):
            if k in rails:
                continue
            poly = np.array([n[1] for n in sp], dtype=float)
            s_a = nearest_s_to_polyline(rails[ia], poly)
            s_b = nearest_s_to_polyline(rails[ib], poly)
            idx_a = piece_index(s_a, cut_s_a)
            idx_b = piece_index(s_b, cut_s_b)
            depth_b = distance_into_piece(s_b, cut_s_b, idx_b, len_b)
            if reversed_b:
                idx_b = n_pieces - 1 - idx_b
            if idx_a == idx_b:
                groups[idx_a].append(sp)
                continue
            straddled += 1
            depth_a = distance_into_piece(s_a, cut_s_a, piece_index(s_a, cut_s_a), rails[ia].length)
            groups[idx_a if depth_a >= depth_b else idx_b].append(sp)

        # A piece is inside a closed shape when the middle of both of its
        # rail pieces is. No shape outline crosses a rail between two cuts,
        # so each rail piece is wholly inside or wholly outside.
        dropped = set()
        half_in = 0
        if polygons:
            bounds_a = [0.0] + cut_s_a + [rails[ia].length]
            bounds_b = [0.0] + cut_s_b + [len_b]
            for i in range(n_pieces):
                j = n_pieces - 1 - i if reversed_b else i
                mid_a = rails[ia].point_tangent((bounds_a[i] + bounds_a[i + 1]) / 2)[0]
                mid_b = rails[ib].point_tangent((bounds_b[j] + bounds_b[j + 1]) / 2)[0]
                in_a = any(bg.point_in_polygon(mid_a, poly) for poly in polygons)
                in_b = any(bg.point_in_polygon(mid_b, poly) for poly in polygons)
                if in_a and in_b:
                    dropped.add(i)
                elif in_a or in_b:
                    half_in += 1

        parent = el.getparent()
        pos = parent.index(el)
        for i in range(n_pieces):
            if i in dropped:
                continue
            piece = el.copy()
            piece.path = make_csp([pieces_a[i], pieces_b[i]] + groups[i]).to_path()
            parent.insert(pos + 1 + i, piece)
            piece.set_random_id(el.get_id())
        parent.remove(el)
        for shape in shapes:
            shape.getparent().remove(shape)

        if half_in:
            inkex.errormsg(
                f"{half_in} piece(s) had one rail inside the shape and the "
                "other outside, so they were kept. Draw the shape all the "
                "way across the column."
            )
        if straddled:
            inkex.errormsg(
                f"{straddled} rung(s) crossed the cut line (one end on each "
                "side). Each went to the side it reaches further into - "
                "check those near the cut."
            )

    def cuts_from_lines(self, lines, shapes):
        columns = [e for e in self.svg.selection.filter(PathElement) if e not in lines]
        if len(columns) != 1:
            raise inkex.AbortExtension(
                "Select exactly one satin column plus the cut line(s) or "
                "closed shape(s). "
                f"Found {len(columns)} satin column(s) in the selection."
            )
        el = columns[0]
        csp, ia, ib, rails = rail_setup(el)

        # Bring each cut line into the satin column's own coordinates.
        to_local = -el.composed_transform()
        line_csp = []
        polygons = []
        for line in lines:
            path = line.path.to_absolute().transform(to_local @ line.composed_transform())
            if line in shapes:
                for sp in path.to_superpath():
                    sp = close_subpath(sp)
                    line_csp.append(sp)
                    polygons.append(shape_polygon(sp))
            else:
                line_csp.extend(path.to_superpath())

        cuts = {}
        for k in (ia, ib):
            hits = [s for s in line_crossings(rails[k], line_csp)
                    if NODE_SNAP < s < rails[k].length - NODE_SNAP]
            if not hits:
                raise inkex.AbortExtension(
                    "The cut line or shape doesn't cross both rails. Draw it "
                    "all the way across the column."
                )
            cuts[k] = [(s, None) for s in hits]

        if len(cuts[ia]) != len(cuts[ib]):
            raise inkex.AbortExtension(
                f"The cut line(s) or shape(s) cross one rail {len(cuts[ia])} "
                f"time(s) and the other {len(cuts[ib])} time(s). Each cut has "
                "to cross each rail exactly once."
            )
        return el, csp, ia, ib, rails, cuts, polygons

    def cuts_from_nodes(self):
        selected = self.options.selected_nodes
        if not selected:
            raise inkex.AbortExtension(
                "Either select the satin column plus a line drawn across it "
                "where you want to cut (or a closed shape drawn over the part "
                "to remove), or, in the Node tool (N), select "
                "the node(s) on the satin column's rails where you want to "
                "cut. Then run this again."
            )

        try:
            parsed = [sel.split(":", 2) for sel in selected]
            ids = {p[0] for p in parsed}
            nodes = [(int(p[1]), int(p[2])) for p in parsed]
        except (ValueError, IndexError):
            raise inkex.AbortExtension(f"Couldn't parse selected node(s) {selected!r}.")

        if len(ids) != 1:
            raise inkex.AbortExtension("All selected nodes must be on the same satin column.")

        el = self.svg.getElementById(ids.pop())
        if el is None or not isinstance(el, PathElement):
            raise inkex.AbortExtension("The selected nodes' object isn't a path.")

        csp, ia, ib, rails = rail_setup(el)

        cuts = {ia: [], ib: []}
        for sp, nl in sorted(set(nodes)):
            if sp not in rails:
                raise inkex.AbortExtension(
                    f"A selected node is on subpath {sp}, which is a rung, not "
                    "a rail. Select nodes on the two rails only."
                )
            rail = rails[sp]
            if nl <= 0 or nl >= len(csp[sp]) - 1:
                raise inkex.AbortExtension(
                    "A selected node is a rail's end node. Cut points have "
                    "to be interior nodes."
                )
            cuts[sp].append((rail.seg_bounds[nl], nl))

        # One rail with no nodes selected: cut it where a perpendicular from
        # each selected node on the other rail crosses it.
        for src, dst in ((ia, ib), (ib, ia)):
            if cuts[src] and not cuts[dst]:
                for s, _ in cuts[src]:
                    p, _ = rails[src].point_tangent(s)
                    d = unit_normal(local_tangent(rails[src], s))
                    res = rails[dst].line_cross(np.array(p), d)
                    if res is None:
                        raise inkex.AbortExtension(
                            "Couldn't find where a line square to the rail at "
                            "a selected node crosses the other rail. Select "
                            "a node on the other rail too."
                        )
                    s_dst = res[1]
                    if s_dst <= NODE_SNAP or s_dst >= rails[dst].length - NODE_SNAP:
                        raise inkex.AbortExtension(
                            "The matching cut on the other rail lands on its "
                            "end. Select a node on the other rail too."
                        )
                    cuts[dst].append((s_dst, None))
                cuts[dst].sort()

        if len(cuts[ia]) != len(cuts[ib]):
            raise inkex.AbortExtension(
                f"One rail has {len(cuts[ia])} selected node(s) and the other "
                f"has {len(cuts[ib])}. Select the same number on each rail "
                "(one per cut), or nodes on only one rail."
            )

        return el, csp, ia, ib, rails, cuts


if __name__ == "__main__":
    CutSatinAtNodes().run()
