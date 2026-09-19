"""
Pure-numpy cubic Bezier geometry helpers: splitting, curve-curve intersection,
polygon flattening, and point-in-polygon testing.

No dependency on inkex, svgpathtools, or shapely - these are not reliably
available in Inkscape's bundled Python environment. Self-contained so it can
be unit-tested with plain `python3 bezier_geom.py`.
"""
import numpy as np

TOL_FLAT = 1e-3      # curve is "flat enough" once control points sit this
                      # close (doc units) to the P0-P3 chord
MAX_DEPTH = 24
DEDUP_T_EPS = 1e-4
DEDUP_PT_EPS = 1e-4


def lerp(a, b, t):
    return a + (b - a) * t


def split_cubic(p0, p1, p2, p3, t):
    """De Casteljau split of a cubic bezier at parameter t. Returns
    (left4, right4), each a tuple of 4 points; left[3] == right[0] == the
    point on the curve at t."""
    p0, p1, p2, p3 = (np.array(p, dtype=float) for p in (p0, p1, p2, p3))
    p01 = lerp(p0, p1, t)
    p12 = lerp(p1, p2, t)
    p23 = lerp(p2, p3, t)
    p012 = lerp(p01, p12, t)
    p123 = lerp(p12, p23, t)
    p0123 = lerp(p012, p123, t)
    to_t = lambda p: (float(p[0]), float(p[1]))
    return (to_t(p0), to_t(p01), to_t(p012), to_t(p0123)), \
           (to_t(p0123), to_t(p123), to_t(p23), to_t(p3))


def _bbox(pts, pad=1e-9):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return (min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad)


def _bbox_overlap(b1, b2):
    return not (b1[2] < b2[0] or b2[2] < b1[0] or b1[3] < b2[1] or b2[3] < b1[1])


def _flatness(p0, p1, p2, p3):
    """Max perpendicular distance of the two control points from the chord
    p0-p3. 0 for a straight/degenerate segment."""
    chord = np.array(p3) - np.array(p0)
    length = np.hypot(*chord)
    if length < 1e-12:
        # Degenerate chord (start==end): fall back to control-point spread.
        return max(np.hypot(*(np.array(p1) - np.array(p0))),
                    np.hypot(*(np.array(p2) - np.array(p0))))
    n = np.array([-chord[1], chord[0]]) / length
    d1 = abs(np.dot(np.array(p1) - np.array(p0), n))
    d2 = abs(np.dot(np.array(p2) - np.array(p0), n))
    return max(d1, d2)


def _line_intersect(a0, a1, b0, b1):
    """Intersect segment a0-a1 with segment b0-b1. Returns (s, u, point) with
    s,u in [0,1] each, or None if parallel/no overlap within [0,1]."""
    a0, a1, b0, b1 = map(np.array, (a0, a1, b0, b1))
    d1 = a1 - a0
    d2 = b1 - b0
    denom = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(denom) < 1e-12:
        return None
    diff = b0 - a0
    s = (diff[0] * d2[1] - diff[1] * d2[0]) / denom
    u = (diff[0] * d1[1] - diff[1] * d1[0]) / denom
    if -1e-9 <= s <= 1 + 1e-9 and -1e-9 <= u <= 1 + 1e-9:
        pt = a0 + d1 * s
        return (min(max(s, 0.0), 1.0), min(max(u, 0.0), 1.0), (pt[0], pt[1]))
    return None


def _bezier_eval_deriv(c, t):
    p0, p1, p2, p3 = (np.array(p, dtype=float) for p in c)
    mt = 1 - t
    pt = (mt**3) * p0 + 3 * (mt**2) * t * p1 + 3 * mt * (t**2) * p2 + (t**3) * p3
    d = 3 * (mt**2) * (p1 - p0) + 6 * mt * t * (p2 - p1) + 3 * (t**2) * (p3 - p2)
    return pt, d


def _refine_t(c, target_pt, t0, iters=12):
    """Newton-Raphson: find t in [0,1] minimizing |eval(c,t) - target_pt|.
    Needed because a chord-fraction (from a straight-line approximation) is
    NOT the same as the true cubic Bezier parameter unless the control
    points are evenly spaced - e.g. inkex represents straight `L` segments
    with both control points collapsed onto the endpoints, which traces a
    straight line but with a non-linear (eased) speed profile."""
    t = t0
    target = np.array(target_pt, dtype=float)
    for _ in range(iters):
        pt, d = _bezier_eval_deriv(c, t)
        denom = np.dot(d, d)
        if denom < 1e-14:
            break
        step = np.dot(pt - target, d) / denom
        new_t = t - step
        if abs(new_t - t) < 1e-13:
            t = new_t
            break
        t = min(max(new_t, 0.0), 1.0)
    return t


def cubic_intersect(c1, c2, depth=0, r1=(0.0, 1.0), r2=(0.0, 1.0)):
    """Recursive subdivision intersection of two cubic beziers.
    c1, c2: ((x,y)*4) control point tuples.
    Returns list of (t1, t2, (x, y)) in the ORIGINAL parameter space of each
    curve (not deduplicated - caller should dedup)."""
    b1, b2 = _bbox(c1), _bbox(c2)
    if not _bbox_overlap(b1, b2):
        return []

    flat1 = _flatness(*c1) < TOL_FLAT
    flat2 = _flatness(*c2) < TOL_FLAT

    if (flat1 and flat2) or depth >= MAX_DEPTH:
        hit = _line_intersect(c1[0], c1[3], c2[0], c2[3])
        if hit is None:
            return []
        s, u, pt = hit
        # s, u are chord fractions, not true Bezier parameters unless control
        # points happen to be evenly spaced - refine against the actual
        # (local, un-shifted) curves before mapping through r1/r2.
        s = _refine_t(c1, pt, s)
        u = _refine_t(c2, pt, u)
        t1 = r1[0] + s * (r1[1] - r1[0])
        t2 = r2[0] + u * (r2[1] - r2[0])
        return [(t1, t2, pt)]

    results = []
    if flat1:
        c1_halves = [(c1, r1)]
    else:
        l, r = split_cubic(*c1, 0.5)
        mid = (r1[0] + r1[1]) / 2
        c1_halves = [(l, (r1[0], mid)), (r, (mid, r1[1]))]

    if flat2:
        c2_halves = [(c2, r2)]
    else:
        l, r = split_cubic(*c2, 0.5)
        mid = (r2[0] + r2[1]) / 2
        c2_halves = [(l, (r2[0], mid)), (r, (mid, r2[1]))]

    for ca, ra in c1_halves:
        for cb, rb in c2_halves:
            results.extend(cubic_intersect(ca, cb, depth + 1, ra, rb))
    return results


def dedup_intersections(hits):
    """hits: list of (t1, t2, pt). Merge near-duplicates (same crossing
    reported by adjacent subdivisions / adjacent shape edges)."""
    out = []
    for t1, t2, pt in sorted(hits, key=lambda h: h[0]):
        merged = False
        for i, (ot1, ot2, opt) in enumerate(out):
            if abs(t1 - ot1) < DEDUP_T_EPS and \
               abs(pt[0] - opt[0]) < DEDUP_PT_EPS and \
               abs(pt[1] - opt[1]) < DEDUP_PT_EPS:
                merged = True
                break
        if not merged:
            out.append((t1, t2, pt))
    return out


def flatten_cubic(p0, p1, p2, p3, n=24):
    """Sample n+1 points along a cubic for polygon approximation."""
    pts = []
    for i in range(n + 1):
        t = i / n
        pts.append(tuple(split_point(p0, p1, p2, p3, t)))
    return pts


def split_point(p0, p1, p2, p3, t):
    p0, p1, p2, p3 = map(np.array, (p0, p1, p2, p3))
    mt = 1 - t
    return (mt**3) * p0 + 3 * (mt**2) * t * p1 + 3 * mt * (t**2) * p2 + (t**3) * p3


def point_in_polygon(pt, poly):
    """Even-odd ray casting. poly: list of (x,y), implicitly closed."""
    x, y = pt
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if ((y1 > y) != (y2 > y)):
            xin = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x < xin:
                inside = not inside
    return inside


if __name__ == "__main__":
    # --- self tests, run with plain `python3 bezier_geom.py` ---
    failures = 0

    def check(name, cond):
        global failures
        status = "OK" if cond else "FAIL"
        if not cond:
            failures += 1
        print(f"[{status}] {name}")

    # 1. Straight line crossing a straight line at a known point.
    c1 = ((0, 0), (0, 0), (10, 10), (10, 10))   # degenerate cubic == line 0,0->10,10
    c2 = ((0, 10), (0, 10), (10, 0), (10, 0))   # line 0,10->10,0
    hits = dedup_intersections(cubic_intersect(c1, c2))
    check("line x line: exactly 1 hit", len(hits) == 1)
    if hits:
        t1, t2, pt = hits[0]
        check(f"line x line: point ~ (5,5), got {pt}",
              abs(pt[0] - 5) < 1e-2 and abs(pt[1] - 5) < 1e-2)
        check(f"line x line: t1 ~ 0.5, got {t1:.4f}", abs(t1 - 0.5) < 1e-2)

    # 2. A horizontal line crossing a real curved cubic (S-curve) twice.
    curve = ((0, 0), (0, 20), (10, -20), (10, 0))
    hline = ((-5, 5), (-5, 5), (15, 5), (15, 5))
    hits2 = dedup_intersections(cubic_intersect(curve, hline))
    check(f"curve x line: >=2 hits, got {len(hits2)}", len(hits2) >= 2)

    # 3. Split-then-recombine sanity: splitting a cubic at t and re-flattening
    #    should reproduce the same point as evaluating the original at t.
    p0, p1, p2, p3 = (0, 0), (5, 15), (15, -15), (20, 0)
    t = 0.37
    left, right = split_cubic(p0, p1, p2, p3, t)
    direct = split_point(p0, p1, p2, p3, t)
    check(f"split_cubic matches split_point at t={t}",
          abs(left[3][0] - direct[0]) < 1e-9 and abs(left[3][1] - direct[1]) < 1e-9)
    check("split_cubic: left[3] == right[0]", left[3] == right[0])

    # 4. Rectangle cut: a straight line crossing a square through two opposite
    #    edges should give exactly 2 intersections.
    square_edges = [
        ((0, 0), (0, 0), (10, 0), (10, 0)),
        ((10, 0), (10, 0), (10, 10), (10, 10)),
        ((10, 10), (10, 10), (0, 10), (0, 10)),
        ((0, 10), (0, 10), (0, 0), (0, 0)),
    ]
    diag_line = ((-5, 5), (-5, 5), (15, 5), (15, 5))
    all_hits = []
    for edge in square_edges:
        all_hits.extend(cubic_intersect(diag_line, edge))
    merged = dedup_intersections(all_hits)
    check(f"line through square: exactly 2 crossings, got {len(merged)}",
          len(merged) == 2)

    # 4b. Regression: inkex represents straight L segments as degenerate
    #     cubics (both control points collapsed onto the endpoints), which
    #     traces a straight line but with a NON-linear (eased) speed
    #     profile. The reported t must be the true Bezier parameter, not the
    #     chord fraction - verify by splitting the original curve at the
    #     reported t and checking the split point lands exactly on the
    #     reported intersection point (this is what the real extension
    #     depends on).
    eased_line = ((0, 10), (0, 10), (50, 10), (50, 10))  # x:0->50 @ y=10, eased
    vertical_at_20 = ((20, -5), (20, -5), (20, 100), (20, 100))
    hit = dedup_intersections(cubic_intersect(eased_line, vertical_at_20))
    check(f"eased line x vertical: 1 hit, got {len(hit)}", len(hit) == 1)
    if hit:
        t1, _t2, pt = hit[0]
        check(f"eased line: intersection point is (20,10), got {pt}",
              abs(pt[0] - 20) < 1e-3 and abs(pt[1] - 10) < 1e-3)
        left, _right = split_cubic(*eased_line, t1)
        check(f"eased line: split_cubic(t1={t1:.4f}) reproduces the point "
              f"(20,10), got {left[3]}",
              abs(left[3][0] - 20) < 1e-3 and abs(left[3][1] - 10) < 1e-3)
        check(f"eased line: t1 is NOT the naive chord fraction 0.4, got {t1:.4f}",
              abs(t1 - 0.4) > 1e-3)

    # 5. point_in_polygon sanity.
    poly = [(0, 0), (10, 0), (10, 10), (0, 10)]
    check("point_in_polygon: center is inside", point_in_polygon((5, 5), poly))
    check("point_in_polygon: outside point is outside", not point_in_polygon((20, 20), poly))

    print()
    if failures:
        print(f"{failures} check(s) FAILED")
        raise SystemExit(1)
    print("all checks passed")
