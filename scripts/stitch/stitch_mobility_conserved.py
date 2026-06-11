#!/usr/bin/env python3
"""Stitch conserved-tier mobility boxplots into a 3-panel supplemental figure.

Layout (A–C, vertical):
  +------------------+
  |        A         |   NHTS Trips     7 x 4 in
  +------------------+
  |        B         |   Air (DB1B)     7 x 4 in
  +------------------+
  |        C         |   SafeGraph Mob. 7 x 4 in
  +------------------+
  Total: 7 x 12 in
"""

import argparse
import os
import subprocess

from svgutils.compose import Figure, SVG, Text
from svg_helpers import Rect

PT_PER_IN  = 72
PANEL_W_IN = 7
PANEL_H_IN = 3
N_PANELS   = 3
TOTAL_W_IN = PANEL_W_IN              # 7
TOTAL_H_IN = PANEL_H_IN * N_PANELS  # 9

LABEL_FONT = "Arial"
LABEL_SIZE = 18
INSET      = 4


def pt(x):
    return x * PT_PER_IN


def add_svg(path, x_off, y_off, native_w_in=PANEL_W_IN):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    scale = pt(PANEL_W_IN) / pt(native_w_in)
    panel = SVG(path)
    panel.scale(scale)
    panel.move(x_off, y_off)
    return panel


def stitch(scenario, out_path):
    dist_dir = f"figs/{scenario}/dist"

    total_w = pt(TOTAL_W_IN)
    total_h = pt(TOTAL_H_IN)

    y_a = 0
    y_b = pt(PANEL_H_IN)
    y_c = pt(PANEL_H_IN * 2)

    elements = [Rect(total_w, total_h, fill="white")]

    elements.append(add_svg(f"{dist_dir}/state_nhts_conserved_boxplot.svg",         0, y_a))
    elements.append(add_svg(f"{dist_dir}/state_air_db1b_conserved_boxplot.svg",      0, y_b))
    elements.append(add_svg(f"{dist_dir}/state_safegraph_conserved_boxplot.svg",     0, y_c))

    for letter, ly in [("A", y_a), ("B", y_b), ("C", y_c)]:
        elements.append(
            Text(letter, INSET, ly + INSET + LABEL_SIZE,
                 size=LABEL_SIZE, weight="bold", font=LABEL_FONT)
        )

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Figure(f"{total_w}pt", f"{total_h}pt", *elements).save(out_path)
    print(f"Saved: {out_path}  ({TOTAL_W_IN:.1f} x {TOTAL_H_IN:.1f} in)")
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
    parser.add_argument("--out", default="manuscript/figures/supp/mobility_conserved.svg")
    parser.add_argument("--pdf", action="store_true")
    parser.add_argument("--pdf-out", default=None)
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
