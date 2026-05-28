#!/usr/bin/env python3
"""Stitch three figure SVGs vertically with A/B/C subfigure labels."""

import argparse
import os
import subprocess

from svgutils.compose import Figure, SVG, Image, Text
from svg_helpers import Rect

PT_PER_IN = 72
PANEL_W_IN = 7
PANEL_H_IN = 4
W = PANEL_W_IN * PT_PER_IN
H = PANEL_H_IN * PT_PER_IN
LABEL_FONT = "Arial"
LABEL_SIZE = 18
LEFT_MARGIN = 24
LETTER_X = 8


def add_panel(path, x_off, y_off):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    ext = os.path.splitext(path)[1].lower()
    if ext == ".svg":
        panel = SVG(path)
        panel.move(x_off, y_off)
        return panel

    if ext in (".jpg", ".jpeg", ".png"):
        panel = Image(W, H, path)
        panel.move(x_off, y_off)
        return panel

    raise ValueError(f"Unsupported file type: {path}")


def stitch(scenario, out_path):
    fig_dir = f"figs/{scenario}"
    panels = [
        ("A", f"{fig_dir}/time/rr_series_poster.svg"),
        ("B", f"{fig_dir}/time/all_pairs_fold_heatmap.svg"),
        ("C", f"{fig_dir}/time_state_networks/network_conserved_combined.jpg"),
    ]

    total_w = W + LEFT_MARGIN
    total_h = H * len(panels)
    elements = [Rect(total_w, total_h, fill="white")]

    for i, (letter, path) in enumerate(panels):
        y_off = i * H
        elements.append(add_panel(path, LEFT_MARGIN, y_off))
        elements.append(
            Text(letter, LETTER_X, y_off + LABEL_SIZE,
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
    parser.add_argument("--out", default="manuscript/figures/geo_time.svg")
    parser.add_argument("--pdf", action="store_true", help="Also export a PDF")
    parser.add_argument("--pdf-out", default=None, help="PDF output path (default: --out with .pdf extension)")
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
