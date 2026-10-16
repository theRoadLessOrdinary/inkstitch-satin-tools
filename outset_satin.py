#!/usr/bin/env python3
"""
Outset / Inset Satin Column

Select one or more satin column objects. The rails are found the way
Ink/Stitch finds them (the subpaths that cross more than two others, else the
two longest), whatever order they are in; every other subpath is a rung. Each rail is pushed away from the other rail by the given
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


def tip_points(a, b, tip):
    """Ends of either rail that touch the other rail: the pointed tips."""
    out = []
    for p, q in ((a, b), (b, a)):
        for e in (p[0], p[-1]):
            if np.min(np.linalg.norm(q - e, axis=1)) <= tip:
                out.append(e)
    return out


class Rail:
    def __init__(self, pts, other_pts):
        self.pts = pts
        tan = tangents(pts)
        normals = np.stack([-tan[:, 1], tan[:, 0]], axis=1)
        # One outward side for the whole rail, by majority vote. Deciding it
        # point by point flips the direction where the rails meet at a tip.
        votes = 0.0
        for p, nrm in zip(pts, normals):
            inward = other_pts[nearest_idx(other_pts, p)] - p
            votes += np.sign(np.dot(nrm, inward))
        self.side = -1.0 if votes > 0 else 1.0
        self.out = self.side * normals
        self.other = other_pts
        self.k = np.ones(len(pts))
        self.dk = np.zeros(len(pts))

    def set_taper(self, d, unit_mm):
        """Per-point share of the distance, 0 to 1, so pointed tips stay
        pointed. Outset fades in over the last FADE_MM from each end that
        meets the other rail. Inset is held to a share of the local width,
        so the rails never cross and a tip keeps its point."""
        pts, other = self.pts, self.other
        width = np.array([np.min(np.linalg.norm(other - p, axis=1)) for p in pts])
        seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        s = np.concatenate([[0.0], np.cumsum(seg)])
        fade = min(FADE_MM / unit_mm, float(s[-1]) / 2.0)
        k = np.ones(len(pts))
        if fade > 0:
            for tp in tip_points(pts, other, TIP_MM / unit_mm):
                x = np.clip(np.linalg.norm(pts - tp, axis=1) / fade, 0.0, 1.0)
                k = np.minimum(k, x * x * (3.0 - 2.0 * x))
        if d < 0:
            # Smoothed along the rail, so a width that jumps (inside a
            # spiral, say) does not make the rail wavy.
            cap = np.minimum(1.0, INSET_SHARE * width / abs(d))
            sigma = SMOOTH_MM / unit_mm
            wgt = np.exp(-((s[:, None] - s[None, :]) ** 2) / (2.0 * sigma * sigma))
            cap = (wgt @ cap) / np.sum(wgt, axis=1)
            k = np.minimum(k, cap)
        self.k = k
        # How fast that share changes per unit of rail length; the handles
        # need it so the rail stays smooth where the distance changes.
        if len(pts) > 2 and s[-1] > 0:
            ds = np.gradient(s)
            ds[ds == 0] = 1.0
            self.dk = np.gradient(k) / ds
        else:
            self.dk = np.zeros(len(pts))

    def tolerance(self, q, tol):
        """Allowed error at q: tol, or less where the column is narrow, so
        the rails cannot touch near a tip."""
        i = nearest_idx(self.pts, q)
        w = float(np.min(np.linalg.norm(self.other - self.pts[i], axis=1)))
        return max(min(tol, 0.15 * w), tol * 0.1)

    def at(self, q):
        """Share of the distance, its rate of change, and the outward
        direction at the rail point nearest q, interpolated between samples
        so they change smoothly rather than in steps."""
        a, b = self.pts[:-1], self.pts[1:]
        ab = b - a
        L2 = np.sum(ab * ab, axis=1)
        L2[L2 == 0] = 1.0
        t = np.clip(np.sum((q - a) * ab, axis=1) / L2, 0.0, 1.0)
        proj = a + t[:, None] * ab
        i = int(np.argmin(np.sum((proj - q) ** 2, axis=1)))
        f = t[i]
        k = self.k[i] * (1 - f) + self.k[i + 1] * f
        dk = self.dk[i] * (1 - f) + self.dk[i + 1] * f
        o = self.out[i] * (1 - f) + self.out[i + 1] * f
        n = np.linalg.norm(o)
        return k, dk, (o / n if n > 1e-12 else self.out[i])

    def target(self, q, d):
        """How far the offset rail should be from the original at q."""
        return abs(d) * self.at(q)[0]

    def dist_to(self, q):
        return float(np.min(np.linalg.norm(self.pts - q, axis=1)))

    def disp(self, q, d):
        k, _, o = self.at(q)
        return d * k * o

    def end_turn(self, p0, p1, p2):
        """Rate the outward normal turns at p0, the end of a cubic whose
        next two control points are p1 and p2. Zero for a missing handle."""
        d1 = 3.0 * (p1 - p0)
        n = np.linalg.norm(d1)
        if n < 1e-9:
            return 0.0
        d2 = 6.0 * (p0 - 2.0 * p1 + p2)
        kappa = (d1[0] * d2[1] - d1[1] * d2[0]) / n ** 3
        return -self.side * kappa

    def unit_out(self, t):
        n = np.linalg.norm(t)
        if n < 1e-9:
            return None
        t = t / n
        return self.side * np.array([-t[1], t[0]])

    def node_disp(self, P, i, d):
        d = d * self.at(P[i][1])[0]
        return self._node_disp(P, i, d)

    def _node_disp(self, P, i, d):
        """Move a node d along the outward normal. At a corner the normals of
        the two sides differ; the node moves to where both offset sides meet
        (a miter), limited to 2 d so a very sharp corner stays close."""
        l, p, r = P[i]
        nl = nr = None
        if i > 0:
            t = p - l
            if np.linalg.norm(t) < 1e-9:
                t = p - P[i - 1][2]
            if np.linalg.norm(t) < 1e-9:
                t = p - P[i - 1][1]
            nl = self.unit_out(t)
        if i < len(P) - 1:
            t = r - p
            if np.linalg.norm(t) < 1e-9:
                t = P[i + 1][0] - p
            if np.linalg.norm(t) < 1e-9:
                t = P[i + 1][1] - p
            nr = self.unit_out(t)
        if nl is None and nr is None:
            return d * self.at(p)[2]
        if nl is None or nr is None:
            return d * (nl if nr is None else nr)
        c = 1.0 + float(np.dot(nl, nr))
        if c < 0.5:
            m = nl + nr
            mn = np.linalg.norm(m)
            return d * 2.0 * (m / mn if mn > 1e-9 else nl)
        return d * (nl + nr) / c

    def offset(self, sp, d):
        """Offset one rail subpath. Each node moves d along the outward
        normal; each handle is scaled by the curvature of its own segment,
        so bends keep their shape and a corner does not blow up."""
        P = [[np.array(c, dtype=float) for c in n] for n in sp]
        out = []
        for i, (l, p, r) in enumerate(P):
            v = self.node_disp(P, i, d)
            k, dk, n_out = self.at(p)
            dl = d * k
            g = d * dk
            kl = kr = 1.0
            if i > 0:
                kl = 1.0 + dl * self.end_turn(p, l, P[i - 1][2])
            if i < len(P) - 1:
                kr = 1.0 + dl * self.end_turn(p, r, P[i + 1][0])
            kl, kr = max(kl, 0.05), max(kr, 0.05)
            np_ = p + v
            lh = (l - p) * kl - g * np.linalg.norm(l - p) * n_out
            rh = (r - p) * kr + g * np.linalg.norm(r - p) * n_out
            out.append([list(np_ + lh), list(np_), list(np_ + rh)])
        return out


def first_crossing(a, b):
    """Where polyline a first crosses polyline b, as (segment index in a,
    point), or None."""
    p, r = a[:-1], np.diff(a, axis=0)
    q, s = b[:-1], np.diff(b, axis=0)
    lo_a, hi_a = np.minimum(a[:-1], a[1:]), np.maximum(a[:-1], a[1:])
    lo_b, hi_b = np.minimum(b[:-1], b[1:]), np.maximum(b[:-1], b[1:])
    for i in range(len(p)):
        near = np.where((hi_b[:, 0] >= lo_a[i, 0]) & (lo_b[:, 0] <= hi_a[i, 0]) &
                        (hi_b[:, 1] >= lo_a[i, 1]) & (lo_b[:, 1] <= hi_a[i, 1]))[0]
        for j in near:
            den = r[i, 0] * s[j, 1] - r[i, 1] * s[j, 0]
            if abs(den) < 1e-12:
                continue
            w = q[j] - p[i]
            t = (w[0] * s[j, 1] - w[1] * s[j, 0]) / den
            u = (w[0] * r[i, 1] - w[1] * r[i, 0]) / den
            if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
                return i, p[i] + t * r[i]
    return None


def segment_hits(p0, p1, b):
    """Indices of the segments of polyline b that segment p0-p1 crosses,
    with the crossing points, in order along b."""
    q, s = b[:-1], np.diff(b, axis=0)
    r = p1 - p0
    den = r[0] * s[:, 1] - r[1] * s[:, 0]
    ok = np.abs(den) > 1e-12
    den = np.where(ok, den, 1.0)
    w = q - p0
    t = (w[:, 0] * s[:, 1] - w[:, 1] * s[:, 0]) / den
    u = (w[:, 0] * r[1] - w[:, 1] * r[0]) / den
    idx = np.where(ok & (t >= 0) & (t <= 1) & (u >= 0) & (u <= 1))[0]
    return idx, p0 + t[idx, None] * r


def _polylines_cross(a, b):
    return first_crossing(a, b) is not None


def rail_indices(samples, lengths):
    """Pick the two rails the way Ink/Stitch does: a rung crosses exactly
    two other subpaths, a rail crosses more. With three subpaths the rails
    cross one. If that does not give exactly two, take the two longest."""
    n = len(samples)
    if n <= 2:
        return list(range(n))
    counts = [sum(_polylines_cross(samples[i], samples[j]) for j in range(n) if j != i)
              for i in range(n)]
    if n == 3:
        cand = [i for i in range(n) if counts[i] == 1]
    else:
        cand = [i for i in range(n) if counts[i] > 2]
    if len(cand) == 2:
        return sorted(cand)
    return sorted(sorted(range(n), key=lambda i: -lengths[i])[:2])


TOLERANCE_MM = 0.05
TIP_MM = 0.3        # rail ends closer than this to the other rail are a pointed tip
FADE_MM = 3.0       # outset grows from 0 at a pointed tip to full over this length
INSET_SHARE = 0.4   # inset never takes more than this share of the local width
SMOOTH_MM = 0.5     # how far along the rail that width limit is averaged
MAX_SPLIT_ROUNDS = 4


def split_segment(sp, i):
    """Split the cubic between nodes i and i+1 at its middle (de Casteljau).
    The shape is unchanged; one node is added."""
    a, b = sp[i], sp[i + 1]
    p0, p1, p2, p3 = (np.array(a[1]), np.array(a[2]), np.array(b[0]), np.array(b[1]))
    m01, m12, m23 = (p0 + p1) / 2, (p1 + p2) / 2, (p2 + p3) / 2
    m012, m123 = (m01 + m12) / 2, (m12 + m23) / 2
    mid = (m012 + m123) / 2
    new_a = [a[0], a[1], list(m01)]
    new_mid = [list(m012), list(mid), list(m123)]
    new_b = [list(m23), b[1], b[2]]
    return sp[:i] + [new_a, new_mid, new_b] + sp[i + 2:]


def split_at(a, b, t):
    """Split the cubic from node a to node b at t. Returns the three nodes
    that replace a and b."""
    p0, p1, p2, p3 = (np.array(a[1]), np.array(a[2]), np.array(b[0]), np.array(b[1]))
    m01, m12, m23 = p0 + (p1 - p0) * t, p1 + (p2 - p1) * t, p2 + (p3 - p2) * t
    m012, m123 = m01 + (m12 - m01) * t, m12 + (m23 - m12) * t
    mid = m012 + (m123 - m012) * t
    return ([a[0], a[1], list(m01)], [list(m012), list(mid), list(m123)],
            [list(m23), b[1], b[2]])


def trim_loops(sp, max_loop):
    """Cut out small loops where an offset rail crosses itself, as happens on
    an inward bend tighter than the offset distance. A loop longer than
    max_loop is left alone (that is the drawing, not the offset)."""
    n = SAMPLES_PER_SEG
    for _ in range(50):
        pts = sample_subpath(sp)
        seg_len = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        cum = np.concatenate([[0.0], np.cumsum(seg_len)])
        found = None
        for i in range(len(pts) - 1):
            far = pts[i + 2:]
            if len(far) < 2:
                break
            idx, xs = segment_hits(pts[i], pts[i + 1], far)
            if len(idx) == 0:
                continue
            j = i + 2 + int(idx[0])
            if cum[j] - cum[i + 1] > max_loop:
                continue
            found = (i, j, xs[0])
            break
        if found is None:
            return sp
        i, j, x = found
        # i and j are polyline segments; map to (bezier segment, t).
        si, sj = i // n, j // n
        ti = (i % n + np.linalg.norm(x - pts[i]) / max(seg_len[i], 1e-12)) / n
        tj = (j % n + np.linalg.norm(x - pts[j]) / max(seg_len[j], 1e-12)) / n
        a1, mid_a, _ = split_at(sp[si], sp[si + 1], min(ti, 1.0))
        _, mid_b, b2 = split_at(sp[sj], sp[sj + 1], min(tj, 1.0))
        joint = [mid_a[0], list(x), mid_b[2]]
        sp = sp[:si] + [a1, joint, b2] + sp[sj + 2:]
    return sp


def cut_rail(sp, i, x, keep_after):
    """Cut a rail at point x on its sampled segment i; keep the part after x
    (or before it)."""
    n = SAMPLES_PER_SEG
    pts = sample_subpath(sp)
    seg = max(float(np.linalg.norm(pts[i + 1] - pts[i])), 1e-12)
    si = min(i // n, len(sp) - 2)
    t = min((i - si * n + np.linalg.norm(x - pts[i]) / seg) / n, 1.0)
    a, mid, b = split_at(sp[si], sp[si + 1], t)
    mid[1] = list(x)
    if keep_after:
        return [[mid[1], mid[1], mid[2]], b] + sp[si + 2:]
    return sp[:si] + [a, [mid[0], mid[1], mid[1]]]


def close_tips(sa, sb, tips, reach):
    """Where the moved rails cross close to a pointed tip, cut both rails at
    the crossing so they meet there, as a point, instead of crossing over."""
    done = []
    for tp in tips:
        if any(np.linalg.norm(tp - q) < reach / 2 for q in done):
            continue
        done.append(tp)
        # Ends of the other tip of the same column limit how far to look.
        others = [np.linalg.norm(tp - q) for q in tips]
        r = min([reach] + [o / 3.0 for o in others if o >= reach / 2])
        pa, pb = sample_subpath(sa), sample_subpath(sb)
        best = None
        for i in range(len(pa) - 1):
            if min(np.linalg.norm(pa[i] - tp), np.linalg.norm(pa[i + 1] - tp)) > r:
                continue
            idx, xs = segment_hits(pa[i], pa[i + 1], pb)
            for j, x in zip(idx, xs):
                dist = float(np.linalg.norm(x - tp))
                if dist <= r and (best is None or dist > best[0]):
                    best = (dist, i, int(j), x)
        if best is None:
            continue
        _, i, j, x = best
        a_tip_at_start = np.linalg.norm(pa[0] - tp) < np.linalg.norm(pa[-1] - tp)
        b_tip_at_start = np.linalg.norm(pb[0] - tp) < np.linalg.norm(pb[-1] - tp)
        sa = cut_rail(sa, i, x, a_tip_at_start)
        sb = cut_rail(sb, j, x, b_tip_at_start)
    return sa, sb


def offset_rail(sp, rail, d, tol):
    """Offset a rail, then split any segment whose offset strays from the
    distance by more than tol and offset again. Only segments that need it
    get an extra node."""
    sp = [list(n) for n in sp]
    for _ in range(MAX_SPLIT_ROUNDS):
        new = rail.offset(sp, d)
        bad = []
        for i in range(len(new) - 1):
            pts = sample_subpath([new[i], new[i + 1]])
            if any(abs(rail.dist_to(q) - rail.target(q, d)) > rail.tolerance(q, tol)
                   for q in pts):
                bad.append(i)
        if not bad:
            break
        for i in reversed(bad):
            sp = split_segment(sp, i)
    else:
        new = rail.offset(sp, d)
    return trim_loops(new, 8.0 * abs(d))


def shift(node, v):
    return [[c[0] + v[0], c[1] + v[1]] for c in node]


def move_rung(sp, pts, rails, rail_pts, d):
    """Each rung node moves with the rail it is nearest to along the rung:
    by that rail's offset at the point where the rung crosses it. A rung
    that does not cross both rails falls back to the nearest rail."""
    hits = []
    for rail, rp in zip(rails, rail_pts):
        hit = first_crossing(pts, rp)
        if hit is None:
            hits = None
            break
        hits.append((hit[0], rail.disp(hit[1], d)))
    new = []
    for k, n in enumerate(sp):
        q = np.array(n[1])
        if hits:
            idx = k * SAMPLES_PER_SEG
            v = min(hits, key=lambda h: abs(h[0] - idx))[1]
        else:
            v = min(rails, key=lambda r: r.dist_to(q)).disp(q, d)
        new.append(shift(n, v))
    return new


def snap_rung(old_sp, new_sp, old_rails, new_rails):
    """Make a straight two-node rung reach the moved rails: each end is put
    on the rung's line where it crosses its moved rail, plus the overshoot it
    had past the original rail. Returns the rung, and False if a rail could
    not be reached."""
    if len(old_sp) != 2:
        return new_sp, True
    a, b = np.array(new_sp[0][1]), np.array(new_sp[1][1])
    u = b - a
    length = np.linalg.norm(u)
    if length < 1e-9:
        return new_sp, True
    u = u / length
    oa, ob = np.array(old_sp[0][1]), np.array(old_sp[1][1])
    ou = (ob - oa) / max(np.linalg.norm(ob - oa), 1e-12)
    far = 4.0 * length + 10.0
    ok = True
    ends = [a, b]
    for which, sign in ((0, -1.0), (1, 1.0)):
        old_end = oa if which == 0 else ob
        new_end = ends[which]
        # The rail this end belongs to: the old rail nearest to the old end.
        ri = min((0, 1), key=lambda r: np.min(np.linalg.norm(old_rails[r] - old_end, axis=1)))
        idx, xs = segment_hits(oa - ou * far, ob + ou * far, old_rails[ri])
        if len(idx) == 0:
            continue
        ox = min(xs, key=lambda x: np.linalg.norm(x - old_end))
        over = float(np.dot(old_end - ox, ou * sign))
        idx, xs = segment_hits(a - u * far, b + u * far, new_rails[ri])
        if len(idx) == 0:
            ok = False
            continue
        nx = min(xs, key=lambda x: np.linalg.norm(x - new_end))
        ends[which] = nx + u * sign * max(over, 0.0)
    out = []
    for n, e in zip(new_sp, ends):
        v = e - np.array(n[1])
        out.append(shift(n, v))
    return out, ok


class OutsetSatin(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--distance", type=float, default=1.0)

    def effect(self):
        sel = [e for e in self.svg.selection if isinstance(e, PathElement)]
        if not sel:
            raise inkex.AbortExtension("Select one or more satin column ladder objects.")

        missed = 0
        for el in sel:
            csp = [list(sp) for sp in el.path.to_absolute().to_superpath()]
            if len(csp) < 2:
                raise inkex.AbortExtension(
                    f"{el.get_id()} needs at least two rail subpaths.")

            m = el.composed_transform()
            scale = abs(m.a * m.d - m.b * m.c) ** 0.5 or 1.0
            d = self.options.distance / scale

            old_csp = [[list(n) for n in sp] for sp in csp]
            samples = [sample_subpath(sp) for sp in csp]
            lengths = [float(np.sum(np.linalg.norm(np.diff(p, axis=0), axis=1)))
                       for p in samples]
            r0, r1 = rail_indices(samples, lengths)
            pts0, pts1 = samples[r0], samples[r1]
            rails = [Rail(pts0, pts1), Rail(pts1, pts0)]
            for rail in rails:
                rail.set_taper(d, scale)

            for ri, si in enumerate((r0, r1)):
                csp[si] = offset_rail(csp[si], rails[ri], d, TOLERANCE_MM / scale)
            csp[r0], csp[r1] = close_tips(csp[r0], csp[r1],
                                          tip_points(pts0, pts1, TIP_MM / scale),
                                          FADE_MM / scale)

            for si in range(len(csp)):
                if si in (r0, r1):
                    continue
                csp[si] = move_rung(csp[si], samples[si], rails, (pts0, pts1), d)

            new_rails = (sample_subpath(csp[r0]), sample_subpath(csp[r1]))
            for si in range(len(csp)):
                if si in (r0, r1):
                    continue
                old = [[list(c) for c in n] for n in old_csp[si]]
                csp[si], ok = snap_rung(old, csp[si], (pts0, pts1), new_rails)
                if not ok:
                    missed += 1

            el.path = make_csp(csp).to_path()

        if missed:
            inkex.errormsg(f"{missed} rung(s) no longer reach both rails, usually next to a "
                           "pointed tip. Check them or delete them.")


if __name__ == "__main__":
    OutsetSatin().run()
