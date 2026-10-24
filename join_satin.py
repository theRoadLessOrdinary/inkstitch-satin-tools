#!/usr/bin/env python3
"""
Join Satin (Ink/Stitch satin column helper)

The reverse of Cut Satin: select two or more satin columns that meet end
to end, where each rail's end node sits on (or very near) a rail end of
the next column, and run this. They become one satin column. Each rail
joins the rail on its own side of the next column, the two end nodes at
each joint become one node (at their midpoint, handles kept), and every
rung is kept.

The columns can be selected in any order and drawn in any direction: the
chain is found from which ends meet. Rails that run opposite ways, or
that swap sides from one column to the next, are lined up first.

The joined column keeps the style and every attribute (Ink/Stitch params
included) of the lowest selected column in the stacking order, and takes
its place; the other columns are removed. Columns in a closed loop are
joined into one column that starts and ends at the same joint.

With nodes selected in the Node tool, only the joints where a selected
node sits are joined; ends that touch elsewhere are left apart. So two
pieces of a loop, joined at one joint, make one column that starts and
ends at the other.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, PathElement

from rung_rail_satin import find_rail_indices, make_csp, reverse_subpath

# Rail ends closer than this count as meeting.
JOIN_TOL_MM = 0.5


def point(node):
    return np.array(node[1], dtype=float)


class Column:
    """One selected column in the first column's coordinates, with rail B
    turned to run the same way as rail A."""

    def __init__(self, el, csp):
        self.el = el
        ia, ib = find_rail_indices(csp)
        a, b = csp[ia], csp[ib]
        same = np.linalg.norm(point(a[0]) - point(b[0])) + np.linalg.norm(point(a[-1]) - point(b[-1]))
        opposite = np.linalg.norm(point(a[0]) - point(b[-1])) + np.linalg.norm(point(a[-1]) - point(b[0]))
        self.b_was_reversed = opposite < same
        self.a = a
        self.b = reverse_subpath(b) if self.b_was_reversed else b
        self.rungs = [sp for i, sp in enumerate(csp) if i not in (ia, ib)]
        self.reversed = False

    def end(self, which):
        """(rail A point, rail B point) at the start (0) or end (1)."""
        k = 0 if which == 0 else -1
        return point(self.a[k]), point(self.b[k])

    def reverse(self):
        self.a, self.b = reverse_subpath(self.a), reverse_subpath(self.b)
        self.reversed = not self.reversed

    def swap(self):
        self.a, self.b = self.b, self.a


def gap(p, q):
    """How far apart two column ends are, with the rails matched either way."""
    straight = max(np.linalg.norm(p[0] - q[0]), np.linalg.norm(p[1] - q[1]))
    crossed = max(np.linalg.norm(p[0] - q[1]), np.linalg.norm(p[1] - q[0]))
    return min(straight, crossed)


def merge_rails(first, second):
    """Join two rails whose last and first nodes meet: one node at the
    midpoint, each side keeping its own handle (moved with the node)."""
    first = [list(map(list, n)) for n in first]
    second = [list(map(list, n)) for n in second]
    p, q = point(first[-1]), point(second[0])
    mid = (p + q) / 2
    in_handle = np.array(first[-1][0]) + (mid - p)
    out_handle = np.array(second[0][2]) + (mid - q)
    joint = [list(in_handle), list(mid), list(out_handle)]
    return first[:-1] + [joint] + second[1:]


class JoinSatin(EffectExtension):
    def effect(self):
        els = [e for e in self.svg.selection.filter(PathElement)
               if len(e.path.to_superpath()) >= 2]
        if len(els) < 2:
            raise inkex.AbortExtension(
                "Select two or more satin columns that meet end to end.")
        doc_order = {e: k for k, e in enumerate(self.svg.iter())}
        els.sort(key=lambda e: doc_order.get(e, 0))
        # Lowest in the stacking order keeps its settings and place.
        base = els[0]
        to_base = -base.composed_transform()
        cols = []
        for el in els:
            path = el.path.to_absolute().transform(to_base @ el.composed_transform())
            csp = [[list(map(list, n)) for n in sp] for sp in path.to_superpath()]
            cols.append(Column(el, csp))

        tol = self.svg.unittouu(f'{JOIN_TOL_MM}mm')
        chain, ends = self.find_chain(cols, tol, self.selected_points(to_base, tol))

        # Walk the chain, turning each column to run onward from the last
        # and putting its rails on the same sides.
        first = chain[0]
        # The first column has to end at its joint with the second, not at
        # its free end.
        if len(chain) > 1 and ends[0][0] == 0:
            first.reverse()
        a, b = first.a, first.b
        rungs = list(first.rungs)
        for k, col in enumerate(chain[1:]):
            if ends[k][1] == 1:
                col.reverse()
            pa, pb = col.end(0)
            if np.linalg.norm(point(a[-1]) - pa) + np.linalg.norm(point(b[-1]) - pb) > \
               np.linalg.norm(point(a[-1]) - pb) + np.linalg.norm(point(b[-1]) - pa):
                col.swap()
            a = merge_rails(a, col.a)
            b = merge_rails(b, col.b)
            rungs += col.rungs

        # Stitch the way the column whose settings are kept did.
        if cols[0].reversed:
            a, b = reverse_subpath(a), reverse_subpath(b)
        if cols[0].b_was_reversed:
            b = reverse_subpath(b)
        base.path = make_csp([a, b] + rungs).to_path()
        # The old node types no longer match the nodes; Inkscape works them out again.
        base.attrib.pop(inkex.addNS('nodetypes', 'sodipodi'), None)
        for col in cols[1:]:
            col.el.getparent().remove(col.el)

    def selected_points(self, to_base, tol):
        """Where the nodes selected in the Node tool are, in the base
        column's coordinates; None when no nodes are selected."""
        points = []
        for sel in self.options.selected_nodes or []:
            try:
                node_id, sp, nl = sel.split(':', 2)
                el = self.svg.getElementById(node_id)
                path = el.path.to_absolute().transform(to_base @ el.composed_transform())
                points.append(point(path.to_superpath()[int(sp)][int(nl)]))
            except (ValueError, IndexError, AttributeError):
                continue
        return points or None

    def find_chain(self, cols, tol, selected=None):
        """The columns in joining order, from which ends meet, and for each
        joint (end of the earlier column, end of the later one). With
        selected points, only ends with a selected node join."""
        n = len(cols)
        links = []
        for i in range(n):
            for j in range(i + 1, n):
                for ei in (0, 1):
                    for ej in (0, 1):
                        d = gap(cols[i].end(ei), cols[j].end(ej))
                        if d <= tol and (selected is None or any(
                                np.linalg.norm(p - q) <= tol for p in selected
                                for q in cols[i].end(ei) + cols[j].end(ej))):
                            links.append((d, i, ei, j, ej))
        links.sort()
        taken = set()
        neighbors = {i: [] for i in range(n)}
        joint = {}
        for _d, i, ei, j, ej in links:
            if (i, ei) in taken or (j, ej) in taken or j in neighbors[i]:
                continue
            taken.update(((i, ei), (j, ej)))
            neighbors[i].append(j)
            neighbors[j].append(i)
            joint[(i, j)], joint[(j, i)] = (ei, ej), (ej, ei)

        # Walk from a column with a free end; in a closed loop, from the
        # lowest one.
        ends = [i for i in range(n) if len(neighbors[i]) < 2]
        start = ends[0] if ends else 0
        order, prev, cur = [start], None, start
        while True:
            nxt = [k for k in neighbors[cur] if k != prev and k not in order]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            order.append(cur)
        if len(order) < n:
            missing = [cols[i].el.get_id() for i in range(n) if i not in order]
            raise inkex.AbortExtension(
                f"These columns don't meet the others end to end: {', '.join(missing)}. "
                f"The end node of each rail has to be within {JOIN_TOL_MM} mm of a rail "
                "end of the next column, and when nodes are selected, a selected node has "
                "to sit at the joint.")
        ends = [joint[(order[k], order[k + 1])] for k in range(len(order) - 1)]
        return [cols[i] for i in order], ends


if __name__ == "__main__":
    JoinSatin().run()
