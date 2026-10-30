#!/usr/bin/env python3
"""
Reroute Rail Along Shape

Select one closed shape plus one or more rails (open paths, or satin
columns). Wherever a rail runs through the shape, the part inside is
replaced by the shape's own edge between the two crossing points, so the
rail now goes around the shape instead of through it. This is what
Path > Difference does to a filled area, applied to an open rail.

- The edge is followed the shorter way around the shape.
- Crossing points are exact (curves are split, not approximated).
- A rail that starts or ends inside the shape is trimmed at the edge.
- A rail that never crosses the shape's edge is left as it is.
- On a satin column (more than 2 subpaths) only the 2 longest subpaths,
  the rails, are rerouted; rungs are left as they are.
- The shape is deleted afterward, as Difference does.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import math

import inkex
from inkex import EffectExtension, ShapeElement, PathElement
from inkex.paths import CubicSuperPath

import bezier_geom as bg

CLOSE_TOL = 1e-6
U_EPS = 1e-6      # rail/ring parameter distance treated as the same spot
MIN_PIECE = 1e-3  # pieces shorter than this (doc units) are grazes, skipped


def make_csp(subpaths):
    obj = CubicSuperPath.__new__(CubicSuperPath)
    list.__init__(obj, subpaths)
    return obj


def subpath_closed(sp, tol=CLOSE_TOL):
    if len(sp) < 2:
        return False
    a, b = sp[0][1], sp[-1][1]
    return abs(a[0] - b[0]) < tol and abs(a[1] - b[1]) < tol


def is_closed(csp):
    return len(csp) > 0 and all(subpath_closed(sp) for sp in csp)


def world_superpath(el):
    return el.path.to_absolute().transform(el.composed_transform()).to_superpath()


def seg(sp, i):
    return (tuple(sp[i][1]), tuple(sp[i][2]), tuple(sp[i + 1][0]), tuple(sp[i + 1][1]))


def sub_cubic(c, ta, tb):
    """The part of cubic c between parameters ta < tb."""
    if tb < 1 - 1e-12:
        c, _ = bg.split_cubic(*c, tb)
        ta = ta / tb if tb > 1e-12 else 0.0
    if ta > 1e-12:
        _, c = bg.split_cubic(*c, ta)
    return c


def piece(sp, u0, u1):
    """Node list ([in, pt, out] triples) of subpath sp from parameter u0 to
    u1 (u = segment index + t), u0 < u1, no wrapping."""
    nseg = len(sp) - 1
    cubics = []
    i = int(math.floor(u0))
    while i < nseg and i < u1 - 1e-12:
        ta = max(u0 - i, 0.0)
        tb = min(u1 - i, 1.0)
        if tb - ta > 1e-12:
            cubics.append(sub_cubic(seg(sp, i), ta, tb))
        i += 1
    return cubics_to_nodes(cubics)


def cubics_to_nodes(cubics):
    if not cubics:
        return []
    nodes = [[list(cubics[0][0]), list(cubics[0][0]), list(cubics[0][1])]]
    for c in cubics:
        nodes[-1][2] = list(c[1])
        nodes.append([list(c[2]), list(c[3]), list(c[3])])
    return nodes


def reverse_nodes(nodes):
    return [[list(n[2]), list(n[1]), list(n[0])] for n in reversed(nodes)]


def ring_arc(ring, v0, v1, forward):
    """Node list along closed ring from parameter v0 to v1, going forward
    (increasing parameter, wrapping) or backward."""
    n = len(ring) - 1
    if not forward:
        return reverse_nodes(ring_arc(ring, v1, v0, True))
    if v1 <= v0:
        v1 += n
    cubics = []
    i = int(math.floor(v0))
    while i < v1 - 1e-12:
        ta = max(v0 - i, 0.0)
        tb = min(v1 - i, 1.0)
        if tb - ta > 1e-12:
            cubics.append(sub_cubic(seg(ring, i % n), ta, tb))
        i += 1
    return cubics_to_nodes(cubics)


def nodes_length(nodes):
    total = 0.0
    for k in range(len(nodes) - 1):
        pts = bg.flatten_cubic(nodes[k][1], nodes[k][2], nodes[k + 1][0], nodes[k + 1][1], n=16)
        for a, b in zip(pts, pts[1:]):
            total += math.hypot(b[0] - a[0], b[1] - a[1])
    return total


def point_at(sp, u):
    nseg = len(sp) - 1
    i = min(int(math.floor(u)), nseg - 1)
    return bg.split_point(*seg(sp, i), u - i)


def subpath_length(sp):
    return nodes_length(sp)


class Shape:
    def __init__(self, csp):
        self.rings = [sp for sp in csp if len(sp) >= 2]
        self.polys = []
        for ring in self.rings:
            pts = []
            for i in range(len(ring) - 1):
                pts.extend(bg.flatten_cubic(*seg(ring, i), n=16)[:-1])
            self.polys.append(pts)

    def inside(self, pt):
        result = False
        for poly in self.polys:
            if bg.point_in_polygon((float(pt[0]), float(pt[1])), poly):
                result = not result
        return result

    def crossings(self, rail):
        """Sorted list of (u on rail, ring index, v on ring) for every place
        the shape's edge crosses the rail."""
        hits = []
        for i in range(len(rail) - 1):
            rs = seg(rail, i)
            for k, ring in enumerate(self.rings):
                for j in range(len(ring) - 1):
                    for t1, t2, pt in bg.cubic_intersect(rs, seg(ring, j)):
                        hits.append((i + t1, k, j + t2, pt))
        hits.sort(key=lambda h: h[0])
        out = []
        for h in hits:
            if out:
                p = out[-1]
                if abs(h[0] - p[0]) < 1e-4 and math.hypot(h[3][0] - p[3][0], h[3][1] - p[3][1]) < 1e-3:
                    continue
            out.append(h)
        return [(u, k, v) for u, k, v, _pt in out]


def join(a, b):
    """Append node list b to a; a's last point and b's first are the same
    spot (a crossing)."""
    if not a:
        return [list(map(list, n)) for n in b]
    if not b:
        return a
    a[-1][2] = list(b[0][2])
    a.extend([list(map(list, n)) for n in b[1:]])
    return a


def reroute_subpath(rail, shape):
    """Returns (list of output subpaths, changed?)."""
    cross = shape.crossings(rail)
    if not cross:
        return [rail], False
    nseg = len(rail) - 1
    stops = [(0.0, None, None)] + cross + [(float(nseg), None, None)]
    out = []
    cur = []
    for a, b in zip(stops, stops[1:]):
        if b[0] - a[0] < U_EPS:
            continue
        part = piece(rail, a[0], b[0])
        if len(part) < 2:
            continue
        mid = point_at(rail, (a[0] + b[0]) / 2.0)
        if not shape.inside(mid) or nodes_length(part) < MIN_PIECE:
            cur = join(cur, part)
            continue
        # Inside piece: swap it for the shape's edge when both ends are
        # crossings on the same ring; otherwise (rail end inside the
        # shape, or ends on different rings) cut it out.
        if a[1] is not None and b[1] is not None and a[1] == b[1]:
            ring = shape.rings[a[1]]
            fwd = ring_arc(ring, a[2], b[2], True)
            back = ring_arc(ring, a[2], b[2], False)
            arc = fwd if nodes_length(fwd) <= nodes_length(back) else back
            if len(arc) >= 2:
                # Pin the arc ends to the rail's own crossing points.
                arc[0][1] = list(part[0][1])
                arc[-1][1] = list(part[-1][1])
                cur = join(cur, arc)
                continue
        if len(cur) >= 2:
            out.append(cur)
        cur = []
    if len(cur) >= 2:
        out.append(cur)
    return out, True


def rail_indices(csp):
    """Which subpaths are rails: all of them, or for a satin column with
    rungs (more than 2 subpaths) the 2 longest."""
    if len(csp) <= 2:
        return set(range(len(csp)))
    lengths = sorted(((subpath_length(sp), i) for i, sp in enumerate(csp)), reverse=True)
    return {lengths[0][1], lengths[1][1]}


class RerouteRail(EffectExtension):
    def effect(self):
        sel = list(self.svg.selection.get(ShapeElement))
        shapes, targets = [], []
        for el in sel:
            try:
                csp = world_superpath(el)
            except (AttributeError, TypeError, ValueError):
                continue
            (shapes if is_closed(csp) else targets).append((el, csp))
        if len(shapes) != 1 or not targets:
            raise inkex.AbortExtension(
                "Select one closed shape and at least one rail (an open path or "
                f"a satin column). Found {len(shapes)} closed shape(s) and "
                f"{len(targets)} open path(s)."
            )
        shape_el, shape_csp = shapes[0]
        shape = Shape(shape_csp)

        changed_any = False
        for el, w_csp in targets:
            rails = rail_indices(w_csp)
            new_subpaths = []
            changed = False
            for i, sp in enumerate(w_csp):
                if i in rails and not subpath_closed(sp) and len(sp) >= 2:
                    parts, ch = reroute_subpath(sp, shape)
                    changed = changed or ch
                    new_subpaths.extend(parts)
                else:
                    new_subpaths.append(sp)
            if not changed:
                continue
            changed_any = True
            self.write_back(el, new_subpaths)

        shape_el.delete()
        if not changed_any:
            inkex.utils.errormsg(
                "The shape's edge doesn't cross any selected rail, so nothing "
                "was rerouted. The shape was still deleted."
            )

    def write_back(self, el, world_subpaths):
        # Back from document coordinates to the element's own coordinates.
        inv = -el.composed_transform()
        local = []
        for sp in world_subpaths:
            nodes = []
            for node in sp:
                nodes.append([list(inv.apply_to_point(p)) for p in node])
            local.append(nodes)
        path = make_csp(local).to_path()
        if isinstance(el, PathElement):
            el.path = path
            return
        new_el = PathElement()
        for key, value in el.attrib.items():
            if key not in ("x", "y", "width", "height", "x1", "y1", "x2", "y2",
                           "points", "cx", "cy", "r", "rx", "ry"):
                new_el.set(key, value)
        new_el.path = path
        el.getparent().replace(el, new_el)


if __name__ == "__main__":
    RerouteRail().run()
