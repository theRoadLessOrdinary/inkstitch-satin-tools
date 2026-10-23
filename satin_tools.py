#!/usr/bin/env python3
"""
Satin Rail Tools (Ink/Stitch satin column toolbox)

Single consolidated Ink/Stitch submenu entry, replacing the 5 separate
entries this extension set used to register (Edit Ladder Rungs, Add Rung
at Node, Cross-Splice Rails, Normalize Stitch Direction, Outset/Inset
Satin Column) with one dialog that has a tab per tool.

This file adds no rung/rail logic of its own - each tab just builds the
matching tool's own EffectExtension subclass, points its .svg at this
extension's already-loaded document (so it sees the same selection) and
its .options at a small stand-in carrying just the option(s) that tool
expects, then calls that tool's unmodified .effect(). The individual
.inx entries and their scripts are left in place and still work
independently; this is purely an additional, consolidated way to reach
the same code.
"""
import os
import sys
from argparse import Namespace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import inkex
from inkex import EffectExtension

from edit_ladder_rungs import EditLadderRungs
from add_rung_at_node import AddRungAtNode
from cross_splice_satin import CrossSpliceSatin
from cut_satin_at_nodes import CutSatinAtNodes
from normalize_stitch_direction import NormalizeStitchDirection
from outset_satin import OutsetSatin


class SatinTools(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--tab", default="rungs")
        pars.add_argument("--rungs_mode", default="add")
        pars.add_argument("--outset_distance", type=float, default=1.0)
        pars.add_argument("--rungs_angle", type=float, default=30.0)

    def effect(self):
        tab = self.options.tab

        if tab == "rungs":
            tool = EditLadderRungs()
            tool.svg = self.svg
            tool.options = Namespace(mode=self.options.rungs_mode, angle=self.options.rungs_angle)

        elif tab == "add_node":
            tool = AddRungAtNode()
            tool.svg = self.svg
            tool.options = Namespace(selected_nodes=self.options.selected_nodes)

        elif tab == "cross_splice":
            tool = CrossSpliceSatin()
            tool.svg = self.svg
            tool.options = Namespace()

        elif tab == "cut":
            tool = CutSatinAtNodes()
            tool.svg = self.svg
            tool.options = Namespace(selected_nodes=self.options.selected_nodes)

        elif tab == "normalize":
            tool = NormalizeStitchDirection()
            tool.svg = self.svg
            tool.options = Namespace()

        elif tab == "outset":
            tool = OutsetSatin()
            tool.svg = self.svg
            tool.options = Namespace(distance=self.options.outset_distance)

        else:
            raise inkex.AbortExtension(f"Unknown tab '{tab}'.")

        tool.effect()


if __name__ == "__main__":
    SatinTools().run()
