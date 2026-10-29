#!/usr/bin/env python3
"""
Cut Satin at Nodes (Ink/Stitch satin column helper)

Either select one or more satin columns (2 rails + rungs, one combined
path each) plus one or more lines drawn across them where you want to cut,
or, in the Node tool (N), select the cut node(s) on the columns' rails.
Every selected column a line or shape crosses is cut; columns it misses
are left alone. Then run this. Each rail is cut at
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

Only satin column(s) selected, with no cut line, shape or nodes: a column
already in separate pieces becomes one satin column per piece (see
separate_sections).

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


def subpath_polyline(sp):
    pts = [np.array(sp[0][1], dtype=float)]
    for j in range(len(sp) - 1):
        pts.extend(np.array(p, dtype=float)
                   for p in bg.flatten_cubic(sp[j][1], sp[j][2], sp[j + 1][0], sp[j + 1][1])[1:])
    return np.array(pts)


def polylines_cross(pa, pb):
    """Whether two polylines meet anywhere other than an end of each. Rail
    pieces left end to end by an earlier cut touch only at their ends, and
    that doesn't join them; a rung crossing a rail does."""
    lo_a, hi_a = pa.min(0), pa.max(0)
    lo_b, hi_b = pb.min(0), pb.max(0)
    if np.any(hi_a < lo_b - NODE_SNAP_DIST) or np.any(hi_b < lo_a - NODE_SNAP_DIST):
        return False
    a0, a1 = pa[:-1, None, :], pa[1:, None, :]
    b0, b1 = pb[None, :-1, :], pb[None, 1:, :]
    da, db = a1 - a0, b1 - b0
    denom = da[..., 0] * db[..., 1] - da[..., 1] * db[..., 0]
    diff = b0 - a0
    with np.errstate(divide='ignore', invalid='ignore'):
        t = (diff[..., 0] * db[..., 1] - diff[..., 1] * db[..., 0]) / denom
        u = (diff[..., 0] * da[..., 1] - diff[..., 1] * da[..., 0]) / denom
    hit = (np.abs(denom) > 1e-12) & (t >= 0) & (t <= 1) & (u >= 0) & (u <= 1)
    if not hit.any():
        return False
    pts = (a0 + t[..., None] * da)[hit]
    ends_a, ends_b = pa[[0, -1]], pb[[0, -1]]
    for p in pts:
        at_a = np.min(np.linalg.norm(ends_a - p, axis=1)) < NODE_SNAP_DIST
        at_b = np.min(np.linalg.norm(ends_b - p, axis=1)) < NODE_SNAP_DIST
        if not (at_a and at_b):
            return True
    return False


def mean_distance(pa, pb):
    """Average distance from points of one polyline to the nearest point of
    the other, both ways: small for the two rails of one piece."""
    d = np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=2)
    return (d.min(1).mean() + d.min(0).mean()) / 2


def pair_pieces(pieces_a, pieces_b):
    """For each piece of rail A, the index of the rail B piece that runs
    closest beside it (closest pairs taken first, each piece used once)."""
    polys_a = [subpath_polyline(sp) for sp in pieces_a]
    polys_b = [subpath_polyline(sp) for sp in pieces_b]
    pairs = sorted((mean_distance(pa, pb) + mean_distance(pb, pa), i, j)
                   for i, pa in enumerate(polys_a) for j, pb in enumerate(polys_b))
    b_for_a = [None] * len(pieces_a)
    used = set()
    for _d, i, j in pairs:
        if b_for_a[i] is None and j not in used:
            b_for_a[i] = j
            used.add(j)
    return b_for_a


def separate_sections(csp):
    """Split one satin column path made of separate pieces into one satin
    column per piece: [rail, rail, rungs...] each. Subpaths that cross are in
    the same piece (a rung crosses both rails of its piece). A piece with no
    rungs is two lone rails, paired by distance."""
    n = len(csp)
    polys = [subpath_polyline(sp) for sp in csp]
    group = list(range(n))

    def find(i):
        while group[i] != i:
            group[i] = group[group[i]]
            i = group[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if find(i) != find(j) and polylines_cross(polys[i], polys[j]):
                group[find(i)] = find(j)
    members = {}
    for i in range(n):
        members.setdefault(find(i), []).append(i)

    sections = [m for m in members.values() if len(m) >= 3]
    small = [m for m in members.values() if len(m) < 3]
    if any(len(m) == 2 for m in small):
        raise inkex.AbortExtension(
            "Part of this column has a rung but only one rail. Check the "
            "column in the Node tool, then run this again.")
    lone = [m[0] for m in small]
    pairs = sorted((mean_distance(polys[i], polys[j]), i, j)
                   for k, i in enumerate(lone) for j in lone[k + 1:])
    used = set()
    for _d, i, j in pairs:
        if i not in used and j not in used:
            used.update((i, j))
            sections.append([i, j])
    if len(used) != len(lone):
        raise inkex.AbortExtension(
            "Part of this column is a single rail with no second rail beside "
            "it, so it can't be a satin column on its own.")

    out = []
    for m in sorted(sections, key=min):
        sub = [csp[i] for i in m]
        ia, ib = find_rail_indices(sub)
        out.append([sub[ia], sub[ib]] + [sp for k, sp in enumerate(sub) if k not in (ia, ib)])
    return out


class CutSatinAtNodes(EffectExtension):
    def effect(self):
        lines = [e for e in self.svg.selection.filter(PathElement)
                 if len(e.path.to_superpath()) == 1]
        lines += list(self.svg.selection.filter(Rectangle, Circle, Ellipse, Polygon))
        if not lines and not self.options.selected_nodes and self.svg.selection:
            self.separate_selected()
            return
        shapes = [e for e in lines if is_closed_cutter(e)]
        # Work out every column's cuts before changing anything, so an error
        # on one column leaves the whole document as it was.
        if lines:
            jobs = self.cuts_from_lines(lines, shapes)
        else:
            jobs = self.cuts_from_nodes()

        half_in = straddled = 0
        for job in jobs:
            h, s = self.cut_column(*job)
            half_in += h
            straddled += s
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

    def cut_column(self, el, csp, ia, ib, rails, cuts, polygons):
        """Split one satin column at its cuts. Returns (pieces half inside a
        shape, rungs that straddled a cut)."""
        pieces_a = split_rail(rails[ia], csp[ia], cuts[ia])
        pieces_b = split_rail(rails[ib], csp[ib], cuts[ib])
        n_pieces = len(pieces_a)
        # Each piece of rail B goes with the piece of rail A it runs closest
        # beside. Comparing where the rails start and end can't tell this on
        # a closed loop (each rail starts and ends at one point), and rails
        # are often drawn anti-parallel.
        b_for_a = pair_pieces(pieces_a, pieces_b)
        a_for_b = {b: a for a, b in enumerate(b_for_a)}
        pieces_b = [pieces_b[b] for b in b_for_a]

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
            idx_b = a_for_b[idx_b]
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
                j = b_for_a[i]
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
            piece.attrib.pop(inkex.addNS('nodetypes', 'sodipodi'), None)
            parent.insert(pos + 1 + i, piece)
            piece.set_random_id(el.get_id())
        parent.remove(el)
        return half_in, straddled

    def separate_selected(self):
        """Only satin column(s) selected, no cut line or nodes: split each
        column that is in separate pieces into one column per piece."""
        columns = list(self.svg.selection.filter(PathElement))
        if not columns:
            raise inkex.AbortExtension(
                "Select a satin column that is in separate pieces, or a satin "
                "column plus a cut line or shape.")
        split = 0
        for el in columns:
            csp = [[list(map(list, n)) for n in sp] for sp in el.path.to_absolute().to_superpath()]
            if len(csp) < 4:
                continue
            sections = separate_sections(csp)
            if len(sections) < 2:
                continue
            parent = el.getparent()
            pos = parent.index(el)
            for i, sub in enumerate(sections):
                piece = el.copy()
                piece.path = make_csp(sub).to_path()
                piece.attrib.pop(inkex.addNS('nodetypes', 'sodipodi'), None)
                parent.insert(pos + 1 + i, piece)
                piece.set_random_id(el.get_id())
            parent.remove(el)
            split += 1
        if not split:
            raise inkex.AbortExtension(
                "The selected column is one piece, so there is nothing to "
                "separate. To cut it, also select a line drawn across it, or "
                "select nodes on its rails in the Node tool.")

    def cuts_from_lines(self, lines, shapes):
        """Every selected satin column the cut lines or shapes cross gets cut.
        Columns they don't touch at all are left alone."""
        columns = [e for e in self.svg.selection.filter(PathElement) if e not in lines]
        if not columns:
            raise inkex.AbortExtension(
                "Select the satin column(s) plus the cut line(s) or closed "
                "shape(s). No satin column was found in the selection."
            )
        jobs = []
        for el in columns:
            name = el.get_id()
            csp, ia, ib, rails = rail_setup(el)

            # Bring each cut line into this satin column's own coordinates.
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
                cuts[k] = [(s, None) for s in line_crossings(rails[k], line_csp)
                           if NODE_SNAP < s < rails[k].length - NODE_SNAP]
            if not cuts[ia] and not cuts[ib]:
                continue
            if not cuts[ia] or not cuts[ib]:
                raise inkex.AbortExtension(
                    f"The cut line or shape crosses only one rail of satin "
                    f"column {name}. Draw it all the way across the column. "
                    "Nothing was changed."
                )
            if len(cuts[ia]) != len(cuts[ib]):
                raise inkex.AbortExtension(
                    f"On satin column {name}, the cut line(s) or shape(s) "
                    f"cross one rail {len(cuts[ia])} time(s) and the other "
                    f"{len(cuts[ib])} time(s). Each cut has to cross each "
                    "rail exactly once. Nothing was changed."
                )
            jobs.append((el, csp, ia, ib, rails, cuts, polygons))

        if not jobs:
            raise inkex.AbortExtension(
                "The cut line or shape doesn't cross any selected satin "
                "column. Draw it all the way across the column(s)."
            )
        return jobs

    def cuts_from_nodes(self):
        """Selected rail nodes, on one or more satin columns. Each column is
        cut at its own nodes."""
        selected = self.options.selected_nodes
        if not selected:
            raise inkex.AbortExtension(
                "Either select the satin column(s) plus a line drawn across "
                "where you want to cut (or a closed shape drawn over the part "
                "to remove), or, in the Node tool (N), select "
                "the node(s) on the satin column's rails where you want to "
                "cut. Then run this again."
            )

        try:
            parsed = [sel.split(":", 2) for sel in selected]
            by_id = {}
            for p in parsed:
                by_id.setdefault(p[0], []).append((int(p[1]), int(p[2])))
        except (ValueError, IndexError):
            raise inkex.AbortExtension(f"Couldn't parse selected node(s) {selected!r}.")

        jobs = []
        for el_id, nodes in by_id.items():
            el = self.svg.getElementById(el_id)
            if el is None or not isinstance(el, PathElement):
                raise inkex.AbortExtension(f"The selected nodes' object {el_id} isn't a path.")
            jobs.append(self.node_cuts(el, nodes))
        return jobs

    def node_cuts(self, el, nodes):
        name = el.get_id()
        csp, ia, ib, rails = rail_setup(el)

        cuts = {ia: [], ib: []}
        for sp, nl in sorted(set(nodes)):
            if sp not in rails:
                raise inkex.AbortExtension(
                    f"A selected node on {name} is on subpath {sp}, which is "
                    "a rung, not a rail. Select nodes on the two rails only."
                )
            rail = rails[sp]
            if nl <= 0 or nl >= len(csp[sp]) - 1:
                raise inkex.AbortExtension(
                    f"A selected node on {name} is a rail's end node. Cut "
                    "points have to be interior nodes."
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
                            f"On {name}, couldn't find where a line square to "
                            "the rail at a selected node crosses the other "
                            "rail. Select a node on the other rail too."
                        )
                    s_dst = res[1]
                    if s_dst <= NODE_SNAP or s_dst >= rails[dst].length - NODE_SNAP:
                        raise inkex.AbortExtension(
                            f"On {name}, the matching cut on the other rail "
                            "lands on its end. Select a node on the other "
                            "rail too."
                        )
                    cuts[dst].append((s_dst, None))
                cuts[dst].sort()

        if len(cuts[ia]) != len(cuts[ib]):
            raise inkex.AbortExtension(
                f"On {name}, one rail has {len(cuts[ia])} selected node(s) "
                f"and the other has {len(cuts[ib])}. Select the same number "
                "on each rail (one per cut), or nodes on only one rail."
            )

        return el, csp, ia, ib, rails, cuts, []


if __name__ == "__main__":
    CutSatinAtNodes().run()
