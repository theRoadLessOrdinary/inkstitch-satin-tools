#!/usr/bin/env python3
"""
Satin Center Line (Ink/Stitch satin column helper)

Select one or more satin columns and run this. Each gets a new path down
the middle of its two rails, with the satin's own style (stroke color,
width and the rest), placed just under the satin in the stacking order.
The satin column is kept as it is. The new path is a plain path, not a
satin column, so Ink/Stitch stitches it as a running stitch. Stroke color
can instead be set to a picked color.

Rungs decide which point on one rail goes with which on the other, as in
Ink/Stitch. The line gets a node wherever either rail or a rung has one,
with smooth handles, plus extra nodes only where a curve needs them to stay
within CENTER_TOL_MM of the true middle.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, PathElement

from rung_rail_satin import Rail, find_rail_indices, normalize_rail_directions, make_csp
from cut_satin_at_nodes import nearest_s_to_polyline, subpath_polyline

# How far the drawn line may stray from the true middle.
CENTER_TOL_MM = 0.05
MAX_SPLITS = 6


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


class Middle:
    """Points and directions down the middle of two rails, by arc length
    on rail A, with rail B's matching position taken from the rungs."""

    def __init__(self, a, b, rungs):
        self.A, self.B = Rail(a), Rail(b)
        anchors = [(0.0, 0.0)]
        for rung in rungs:
            poly = subpath_polyline(rung)
            anchors.append((nearest_s_to_polyline(self.A, poly),
                            nearest_s_to_polyline(self.B, poly)))
        anchors.append((self.A.length, self.B.length))
        anchors.sort()
        # Keep only rungs that go forward along both rails.
        kept = [anchors[0]]
        for sa, sb in anchors[1:]:
            if sa > kept[-1][0] + 1e-9 and sb > kept[-1][1] + 1e-9:
                kept.append((sa, sb))
        if kept[-1] != anchors[-1]:
            kept[-1] = anchors[-1]
        self.anchor_a = np.array([k[0] for k in kept])
        self.anchor_b = np.array([k[1] for k in kept])

    def b_at(self, sa):
        return float(np.interp(sa, self.anchor_a, self.anchor_b))

    def a_at(self, sb):
        return float(np.interp(sb, self.anchor_b, self.anchor_a))

    def at(self, sa):
        pa, ta = self.A.point_tangent(sa)
        pb, tb = self.B.point_tangent(self.b_at(sa))
        mid = (np.array(pa, dtype=float) + np.array(pb, dtype=float)) / 2
        return mid, unit(unit(np.array(ta, dtype=float)) + unit(np.array(tb, dtype=float)))

    def stations(self):
        """Arc positions on rail A where the line gets a node: every node of
        either rail and every rung."""
        s = set(self.anchor_a.tolist())
        s.update(self.A.seg_bounds)
        s.update(self.a_at(sb) for sb in self.B.seg_bounds)
        s = sorted(min(max(v, 0.0), self.A.length) for v in s)
        out = [s[0]]
        for v in s[1:]:
            if v - out[-1] > 1e-6 * max(self.A.length, 1.0):
                out.append(v)
        return out


def bezier_point(p0, p1, p2, p3, t):
    mt = 1 - t
    return mt ** 3 * p0 + 3 * mt ** 2 * t * p1 + 3 * mt * t ** 2 * p2 + t ** 3 * p3


def center_line(middle, tol):
    stations = middle.stations()
    nodes = [middle.at(s) for s in stations]

    def segment(s0, n0, s1, n1, depth):
        (p0, d0), (p3, d3) = n0, n1
        k = np.linalg.norm(p3 - p0) / 3
        c1, c2 = p0 + d0 * k, p3 - d3 * k
        sm = (s0 + s1) / 2
        nm = middle.at(sm)
        if depth < MAX_SPLITS and np.linalg.norm(bezier_point(p0, c1, c2, p3, 0.5) - nm[0]) > tol:
            return segment(s0, n0, sm, nm, depth + 1) + segment(sm, nm, s1, n1, depth + 1)
        return [(c1, c2, p3)]

    pieces = []
    for i in range(len(stations) - 1):
        pieces += segment(stations[i], nodes[i], stations[i + 1], nodes[i + 1], 0)

    start = nodes[0][0]
    sp = [[list(start), list(start), list(start)]]
    for c1, c2, p3 in pieces:
        sp[-1][2] = list(c1)
        sp.append([list(c2), list(p3), list(p3)])
    return sp


class SatinCenterLine(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--color_source", default="own", choices=["own", "pick"])
        pars.add_argument("--stroke_color", type=inkex.Color, default=inkex.Color("#000000"))

    def effect(self):
        els = [e for e in self.svg.selection.filter(PathElement)
               if len(e.path.to_superpath()) >= 2]
        if not els:
            raise inkex.AbortExtension("Select one or more satin columns.")
        tol = self.svg.unittouu(f'{CENTER_TOL_MM}mm')
        for el in els:
            csp = [[list(map(list, n)) for n in sp]
                   for sp in el.path.to_absolute().to_superpath()]
            ia, ib = find_rail_indices(csp)
            a, b = normalize_rail_directions(csp[ia], csp[ib])
            rungs = [sp for i, sp in enumerate(csp) if i not in (ia, ib)]
            # Tolerance is on the page; the path's own units may be scaled.
            t = el.composed_transform()
            scale = np.sqrt(abs(t.a * t.d - t.b * t.c)) or 1.0
            sp = center_line(Middle(a, b, rungs), tol / scale)

            line = PathElement()
            line.path = make_csp([sp]).to_path()
            if el.get('style'):
                line.set('style', el.get('style'))
            if self.options.color_source == 'pick':
                # The Tool Panel sends the color as a "#rrggbb" string.
                line.style['stroke'] = str(inkex.Color(self.options.stroke_color).to_rgb())
            if el.get('transform'):
                line.set('transform', el.get('transform'))
            line.set('id', self.svg.get_unique_id('center'))
            line.label = f"{el.label or el.get_id()} center line"
            el.addprevious(line)


if __name__ == "__main__":
    SatinCenterLine().run()
