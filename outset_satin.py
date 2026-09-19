#!/usr/bin/env python3
"""
Outset / Inset Satin Column

Select one or more already-built ladder objects (2 rail subpaths followed by
rung subpaths). Each rail is pushed away from the other rail by the given
distance (negative distance insets, i.e. thinner). Every rung endpoint moves
by the same displacement as the rail point nearest to it, so rungs stay
attached to the rails and keep their overshoot.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, PathElement

from rung_rail_satin import make_csp

SAMPLES_PER_SEG = 24


def _cubic(p0, p1, p2, p3, t):
    u = 1.0 - t
    return (u ** 3)[:, None] * p0 + (3 * u * u * t)[:, None] * p1 + \
        (3 * u * t * t)[:, None] * p2 + (t ** 3)[:, None] * p3


def sample_subpath(sp):
    ts = np.linspace(0.0, 1.0, SAMPLES_PER_SEG, endpoint=False)
    pts = []
    for a, b in zip(sp[:-1], sp[1:]):
        pts.append(_cubic(np.array(a[1]), np.array(a[2]), np.array(b[0]), np.array(b[1]), ts))
    pts.append(np.array([sp[-1][1]]))
    return np.vstack(pts)


def tangents(pts):
    t = np.gradient(pts, axis=0)
    n = np.linalg.norm(t, axis=1)
    n[n == 0] = 1.0
    return t / n[:, None]


def nearest_idx(pts, q):
    return int(np.argmin(np.sum((pts - q) ** 2, axis=1)))


class Rail:
    def __init__(self, pts, other_pts):
        self.pts = pts
        tan = tangents(pts)
        normals = np.stack([-tan[:, 1], tan[:, 0]], axis=1)
        out = np.zeros_like(normals)
        for i, (p, nrm) in enumerate(zip(pts, normals)):
            inward = other_pts[nearest_idx(other_pts, p)] - p
            out[i] = -nrm if np.dot(nrm, inward) > 0 else nrm
        self.out = out

    def dist_to(self, q):
        return float(np.min(np.linalg.norm(self.pts - q, axis=1)))

    def disp(self, q, d):
        return d * self.out[nearest_idx(self.pts, q)]


def shift(node, v):
    return [[c[0] + v[0], c[1] + v[1]] for c in node]


class OutsetSatin(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--distance", type=float, default=1.0)

    def effect(self):
        sel = [e for e in self.svg.selection if isinstance(e, PathElement)]
        if not sel:
            raise inkex.AbortExtension("Select one or more satin column ladder objects.")

        for el in sel:
            csp = [list(sp) for sp in el.path.to_absolute().to_superpath()]
            if len(csp) < 2:
                raise inkex.AbortExtension(
                    f"{el.get_id()} needs at least two rail subpaths.")

            m = el.transform
            scale = abs(m.a * m.d - m.b * m.c) ** 0.5 or 1.0
            d = self.options.distance / scale

            pts0, pts1 = sample_subpath(csp[0]), sample_subpath(csp[1])
            rails = [Rail(pts0, pts1), Rail(pts1, pts0)]

            for ri in (0, 1):
                csp[ri] = [shift(n, rails[ri].disp(np.array(n[1]), d)) for n in csp[ri]]

            for si in range(2, len(csp)):
                new = []
                for n in csp[si]:
                    q = np.array(n[1])
                    rail = min(rails, key=lambda r: r.dist_to(q))
                    new.append(shift(n, rail.disp(q, d)))
                csp[si] = new

            el.path = make_csp(csp).to_path()


if __name__ == "__main__":
    OutsetSatin().run()
