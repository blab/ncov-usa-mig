#!/usr/bin/env python3
"""Stitch the two age-time trace figures into a single stacked SVG/PDF.

Layout (subfigure letters A-B), both panels 9 x 4 in:
  +-----------------------------------+
  | A | vaccine-rollout elderly nRR   |
  +-----------------------------------+
  | B | school-age nRR traces         |
  +-----------------------------------+
"""

import argparse
import os

from svgutils.compose import Figure, SVG, Text
from svg_helpers import Rect


PT_PER_IN = 72
PANEL_W_IN = 9
PANEL_H_IN = 4
W = PANEL_W_IN * PT_PER_IN

LABEL_FONT = "Arial"
LABEL_SIZE = 18
LEFT_MARGIN = 24
LETTER_X = 8


def in_pt(x):
    return x * PT_PER_IN


def add_panel(path, x_off, y_off):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    panel = SVG(path)
    panel.move(x_off, y_off)
    return panel


def stitch(scenario, out_path):
    fig_dir = f"figs/{scenario}/age_time"
    vaccine = f"{fig_dir}/all_countries/vaccine_period_nRR_fixed_elderly_ci.svg"
    school = f"{fig_dir}/school/school_nRR_trace_simplified.svg"

    total_w = W + LEFT_MARGIN
    total_h = in_pt(2 * PANEL_H_IN)

    elements = [
        Rect(total_w, total_h, fill="white"),
        add_panel(vaccine, LEFT_MARGIN, 0),
        add_panel(school, LEFT_MARGIN, in_pt(PANEL_H_IN)),
        Text("A", LETTER_X, LABEL_SIZE, size=LABEL_SIZE, weight="bold", font=LABEL_FONT),
        Text("B", LETTER_X, in_pt(PANEL_H_IN) + LABEL_SIZE,
             size=LABEL_SIZE, weight="bold", font=LABEL_FONT),
    ]

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Figure(f"{total_w}pt", f"{total_h}pt", *elements).save(out_path)
    print(f"Stitched figure saved to: {out_path} "
          f"({total_w/PT_PER_IN:.2f}in x {total_h/PT_PER_IN:.2f}in)")
    return out_path


def export_pdf(svg_path, pdf_path=None):
    import subprocess
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
    parser.add_argument("--out", default="manuscript/figures/age_time.svg")
    parser.add_argument("--pdf", action="store_true", help="Also export a PDF")
    parser.add_argument("--pdf-out", default=None)
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
