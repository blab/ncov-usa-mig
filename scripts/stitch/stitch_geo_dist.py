#!/usr/bin/env python3
"""Stitch geographic distance figure panels into a single multi-panel SVG.

Layout (subfigure letters A-B):
  +------------------------------------------+
  | A | state_distance_poster      (10x2.5in) |
  +------------------------------------------+
  | B | pred_rr_curves_combined    (10x4in)   |
  +------------------------------------------+
"""

import argparse
import os
import subprocess

from svgutils.compose import Figure, SVG, Text
from svg_helpers import Rect

PT_PER_IN  = 72
TOTAL_W_IN = 10
LEFT_MARGIN = 24
LETTER_X   = 8
LABEL_FONT = "Arial"
LABEL_SIZE = 18

# Native dimensions of each source figure (inches)
DIST_W_IN, DIST_H_IN  = 10, 2.5
COMB_W_IN, COMB_H_IN  = 10, 4

TOTAL_H_IN = DIST_H_IN + COMB_H_IN   # 6.5


def in_pt(x):
    return x * PT_PER_IN


def add_svg(path, x_off, y_off, scale=None):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    panel = SVG(path)
    if scale is not None:
        panel.scale(scale)
    panel.move(x_off, y_off)
    return panel


def stitch(scenario, out_path):
    fig_dir = f"figs/{scenario}/dist"

    total_w = in_pt(TOTAL_W_IN) + LEFT_MARGIN
    total_h = in_pt(TOTAL_H_IN)

    elements = [Rect(total_w, total_h, fill="white")]

    # Row A: state_distance_poster — native width matches TOTAL_W_IN, no scaling
    elements.append(add_svg(
        f"{fig_dir}/state_distance_poster.svg",
        LEFT_MARGIN, 0
    ))

    # Row B: pred_rr_curves_combined — native width matches TOTAL_W_IN, no scaling
    elements.append(add_svg(
        f"{fig_dir}/pred_rr_curves_combined.svg",
        LEFT_MARGIN, in_pt(DIST_H_IN)
    ))

    # Subfigure letters
    INSET  = 4
    row2_y = in_pt(DIST_H_IN)
    letters = [
        ("A", LETTER_X, INSET),
        ("B", LETTER_X, row2_y + INSET),
    ]
    for letter, lx, ly in letters:
        elements.append(
            Text(letter, lx, ly + LABEL_SIZE,
                 size=LABEL_SIZE, weight="bold", font=LABEL_FONT)
        )

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Figure(f"{total_w}pt", f"{total_h}pt", *elements).save(out_path)
    print(f"Stitched figure saved to: {out_path} "
          f"({total_w/PT_PER_IN:.2f}in x {total_h/PT_PER_IN:.2f}in)")
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
    parser.add_argument("--out", default="manuscript/figures/geo_dist.svg")
    parser.add_argument("--pdf", action="store_true", help="Also export a PDF")
    parser.add_argument("--pdf-out", default=None, help="PDF output path (default: --out with .pdf extension)")
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
