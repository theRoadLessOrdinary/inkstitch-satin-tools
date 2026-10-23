#!/usr/bin/env python3
"""Subpath numbers for Ladder Rungs.

Number: on each selected path, writes one mapsedge:sbNNNN attribute per
subpath (numbered from 1 in the order they appear in d), whose value is
that subpath's own coordinates, and puts a tiny number beside each subpath
in a "Subpath numbers" layer. Each number carries mapsedge:of (the path id)
and mapsedge:sb, so the number, the attribute and the subpath point at each
other. The layer has Ink/Stitch's "ignore layer" command, so the numbers
never stitch. Numbering again replaces the old numbers.

Remove numbers: deletes the attributes and the numbers for the selected
paths, or for every path when nothing is selected.

The numbers are a snapshot: after the path is edited they go out of date
until it is numbered again.
"""
import os

import inkex
from inkex import PathElement, TextElement, Group, Use
from inkex.paths import Path, Move, move
from lxml import etree

NS = "https://mapsedgecreative.com/ns/inkscape"
etree.register_namespace("mapsedge", NS)
inkex.NSS["mapsedge"] = NS

LAYER_ID = "mapsedge_subpath_numbers"
LAYER_LABEL = "Subpath numbers"
LABEL_MM = 0.8
LABEL_STYLE = "font-family:sans-serif;font-weight:bold;fill:#d000d0;stroke:none;text-anchor:middle"
IGNORE_SYMBOL = "inkstitch_ignore_layer"
INKSTITCH_SYMBOLS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "inkstitch", "symbols", "inkstitch.svg")


def q(name):
    return f"{{{NS}}}{name}"


def subpaths_of(el):
    """Each subpath of el's d, in order, as its own absolute path string."""
    parts = []
    for seg in el.path.to_absolute():
        if isinstance(seg, (Move, move)) or not parts:
            parts.append([])
        parts[-1].append(seg)
    return [str(Path(p)) for p in parts]


def old_numbers(el):
    return {k: v for k, v in el.attrib.items() if k.startswith(f"{{{NS}}}sb")}


def labels_for(svg, pid=None):
    layer = svg.getElementById(LAYER_ID)
    if layer is None:
        return []
    return [t for t in layer if t.get(q("of")) is not None and (pid is None or t.get(q("of")) == pid)]


def number_layer(svg):
    layer = svg.getElementById(LAYER_ID)
    if layer is not None:
        return layer
    layer = Group.new(LAYER_LABEL)
    layer.set("id", LAYER_ID)
    layer.set("inkscape:groupmode", "layer")
    svg.append(layer)  # on top, so the numbers show over the work
    if svg.getElementById(IGNORE_SYMBOL) is None:
        try:
            sym = inkex.load_svg(INKSTITCH_SYMBOLS).getroot().getElementById(IGNORE_SYMBOL)
        except (OSError, etree.XMLSyntaxError):
            sym = None
        if sym is not None:
            svg.defs.append(sym)
    use = Use()
    use.set("xlink:href", f"#{IGNORE_SYMBOL}")
    use.set("x", "0")
    use.set("y", "0")
    use.set("transform", "scale(0.2)")
    layer.append(use)
    return layer


def midpoint(sub_d):
    """A point halfway along one subpath."""
    csp = Path(sub_d).to_superpath()[0]
    lens = [inkex.bezier.cspseglength(a, b) for a, b in zip(csp, csp[1:])]
    if not lens or sum(lens) <= 0:
        return csp[0][1]
    half = sum(lens) / 2.0
    for (a, b), L in zip(zip(csp, csp[1:]), lens):
        if half <= L and L > 0:
            t = inkex.bezier.beziertatlength((a[1], a[2], b[0], b[1]), half / L)
            return inkex.bezier.bezierpointatt((a[1], a[2], b[0], b[1]), t)
        half -= L
    return csp[-1][1]


def number(svg, els):
    """Number the subpaths of els. Returns notes for the user."""
    notes = []
    size = svg.viewport_to_unit(f"{LABEL_MM}mm")
    layer = number_layer(svg)
    for el in els:
        pid = el.get_id()
        subs = subpaths_of(el)
        old = old_numbers(el)
        if old:
            moved = sum(1 for v in old.values() if v not in subs)
            if moved:
                notes.append(f"{el.label or pid}: {moved} of the old numbers no longer matched the path; renumbered.")
        for k in old:
            del el.attrib[k]
        for t in labels_for(svg, pid):
            t.getparent().remove(t)
        to_doc = el.composed_transform()
        for i, sub in enumerate(subs, 1):
            sb = f"{i:04d}"
            el.set(q(f"sb{sb}"), sub)
            x, y = to_doc.apply_to_point(midpoint(sub))
            t = TextElement()
            t.text = str(i)
            t.set("x", f"{x:.4f}")
            t.set("y", f"{y + size * 0.35:.4f}")
            t.set("style", f"{LABEL_STYLE};font-size:{size:.4f}px")
            t.set(q("of"), pid)
            t.set(q("sb"), sb)
            layer.append(t)
    return notes


def unnumber(svg, els):
    """Remove the numbers from els, or from every path when els is empty."""
    if not els:
        els = [e for e in svg.descendants() if isinstance(e, PathElement) and old_numbers(e)]
        for t in labels_for(svg):
            t.getparent().remove(t)
    for el in els:
        for k in old_numbers(el):
            del el.attrib[k]
        for t in labels_for(svg, el.get_id()):
            t.getparent().remove(t)
    layer = svg.getElementById(LAYER_ID)
    if layer is not None and not labels_for(svg):
        layer.getparent().remove(layer)
