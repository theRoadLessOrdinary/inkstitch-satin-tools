#!/usr/bin/env python3
"""
Normalize Stitch Direction (Ink/Stitch satin column helper)

Select two or more already-built ladder (2-rail satin column) objects, in
whatever order they'll actually be stitched (their document/z-order).

The first object's rails get oriented bottom-to-top (in document Y - SVG Y
increases downward, so "bottom" is the larger-Y end). Every following
object then gets oriented (both its rails reversed together, if needed) so
that its start lands as close as possible to wherever the previous object's
rail ended - so the whole selection stitches through in one continuous
direction, beginning to end, instead of jumping to a far point and
doubling back.

Only the two rail subpaths of each object are ever reversed; any rung
subpaths already present are left untouched (their own geometry doesn't
depend on direction).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, PathElement

from rung_rail_satin import reverse_subpath, make_csp


class NormalizeStitchDirection(EffectExtension):
    def effect(self):
        sel = [e for e in self.svg.selection if isinstance(e, PathElement)]
        if len(sel) < 2:
            raise inkex.AbortExtension(
                "Select two or more already-built ladder objects, in stitch order."
            )

        doc_order = {el: i for i, el in enumerate(self.svg.iter())}
        sel.sort(key=lambda e: doc_order.get(e, 0))

        prev_end = None
        for i, el in enumerate(sel):
            csp = el.path.to_absolute().to_superpath()
            if len(csp) < 2:
                raise inkex.AbortExtension(
                    f"{el.get_id()} does not have two rail subpaths - build its rungs first."
                )
            rail_start = np.array(csp[0][0][1])
            rail_end = np.array(csp[0][-1][1])

            if i == 0:
                # bottom-to-top: start should be the larger-Y (lower) point.
                reverse = rail_start[1] < rail_end[1]
            else:
                reverse = (np.linalg.norm(rail_end - prev_end) <
                           np.linalg.norm(rail_start - prev_end))

            if reverse:
                csp[0] = reverse_subpath(csp[0])
                csp[1] = reverse_subpath(csp[1])
                rail_start, rail_end = rail_end, rail_start
                el.path = make_csp(list(csp)).to_path()

            prev_end = rail_end


if __name__ == "__main__":
    NormalizeStitchDirection().run()
