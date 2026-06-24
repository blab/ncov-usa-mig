#!/usr/bin/env python3
"""Stitch descriptive figure panels into a single multi-panel SVG.

Layout (subfigure letters a-d), dimensions from original Inkscape assembly:
  +----------------------+------------------------------+
  | a  seq_map           | b  natl_effort               |
  |    (left, full ht)   |    (top-right, 6x3in → 2:1)  |
  |                      +-----------+------------------+
  |                      | c age_hist| d  age_region    |
  |                      |  (3x4in)  |    (5x4in)       |
  +----------------------+-----------+------------------+

Source SVG dimensions (set in descriptive_figures.R):
  seq_map:     6x7in   age_hist:   3x4in
  natl_effort: 6x3in   age_region: 5x4in
"""

import argparse
import os
import subprocess

from svgutils.compose import Figure, Image, SVG, Text
from svg_helpers import Rect

PT_PER_IN = 72
TARGET_W_IN = 7

# Original Inkscape canvas (mm) — proportions preserved, scale to TARGET_W_IN
ORIG_W_MM = 91.032816
ORIG_H_MM = 50.834633

# pt per mm at target width
MM_TO_PT = (TARGET_W_IN * PT_PER_IN) / ORIG_W_MM

# Panel slots (mm from original Inkscape layout) + native ggsave width (in)
PANELS = {
    "seq_map":     dict(x=-2.1166668, y=0,           w=43.542858, h=50.799999, native_w_in=6),
    "natl_effort": dict(x=40.110756,  y=0.094874,    w=50.799999, h=25.4,      native_w_in=6),
    "age_hist":    dict(x=39.077778,  y=25.434631,   w=19.049999, h=25.4,      native_w_in=3),
    "age_region":  dict(x=59.282818,  y=25.036949,   w=31.75,     h=25.4,      native_w_in=5),
}

# Subfigure labels (mm, from original Inkscape SVG)
LABELS = [
    ("A", 1.9594771,  3.5108116),
    ("B", 37.673016,  3.5108116),
    ("C", 37.741917,  27.988649),
    ("D", 58.678082,  27.988649),
]
LABEL_FONT = "Arial"
LABEL_SIZE = 10


def mm(v):
    """Convert mm (original layout) to pt at target width."""
    return v * MM_TO_PT


def add_svg(path, slot):
    """Load SVG, scale uniformly to fit slot width, position at slot origin."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    target_w_pt = mm(slot["w"])
    native_w_pt = slot["native_w_in"] * PT_PER_IN
    scale = target_w_pt / native_w_pt
    panel = SVG(path)
    panel.scale(scale)
    panel.move(mm(slot["x"]), mm(slot["y"]))
    return panel


def stitch(scenario, out_path):
    fig_dir = f"figs/{scenario}/desc"

    total_w = mm(ORIG_W_MM)
    total_h = mm(ORIG_H_MM)

    elements = [Rect(total_w, total_h, fill="white")]

    for name, slot in PANELS.items():
        if name == "seq_map":
            w_pt = mm(slot["w"])
            h_pt = mm(slot["h"])
            panel = Image(w_pt, h_pt, os.path.join(fig_dir, f"{name}.png"))
            panel.move(mm(slot["x"]), mm(slot["y"]))
            elements.append(panel)
        else:
            svg_path = os.path.join(fig_dir, f"{name}.svg")
            elements.append(add_svg(svg_path, slot))

    for letter, lx_mm, ly_mm in LABELS:
        elements.append(
            Text(letter, mm(lx_mm), mm(ly_mm),
                 size=LABEL_SIZE, weight="bold", font=LABEL_FONT)
        )

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Figure(f"{total_w}pt", f"{total_h}pt", *elements).save(out_path)
    print(f"Stitched figure saved to: {out_path} "
          f"({total_w/PT_PER_IN:.3f}in x {total_h/PT_PER_IN:.3f}in)")
    return out_path


def export_pdf(svg_path, pdf_path=None):
    if pdf_path is None:
        pdf_path = os.path.splitext(svg_path)[0] + ".pdf"
    subprocess.run(
        ["inkscape", "--export-type=pdf", f"--export-filename={pdf_path}", svg_path],
        check=True,
    )
    print(f"PDF exported to: {pdf_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario", default="CAM_1000")
    parser.add_argument("--out", default="manuscript/figures/desc.svg")
    parser.add_argument("--pdf", action="store_true", help="Also export a PDF")
    parser.add_argument("--pdf-out", default=None)
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
