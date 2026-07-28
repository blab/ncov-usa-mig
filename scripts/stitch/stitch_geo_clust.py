#!/usr/bin/env python3
"""Stitch geographic clustering figure panels into a single multi-panel SVG.

Layout:
  +------------------+
  |        A         |   8.000 x 2.133 in  patched_maps.jpg              (raster)
  +------------------+
  |        B         |   8.000 x 6.667 in  state_heatmap_clustered.svg   (vector)
  +------------------+

(The region map + PCoA panels that used to sit below this as C/D now live
in their own figure -- see stitch_geo_pc.py -- since this one was getting
too tall.)
"""

import argparse
import os
import subprocess

from svgutils.compose import Figure, Image, SVG, Text
from svg_helpers import Rect

PT_PER_IN = 72

# Total figure width
TARGET_W_IN = 8.0

# Panel heights (both rows are full width)
H_A_IN = 2.133   # patched_maps scaled to full width
H_B_IN = 6.667   # heatmap scaled to full width

# Native ggsave width for the heatmap svg (used to compute its scale factor)
NATIVE_W_B_IN = 12   # state_heatmap_clustered.svg

LABEL_FONT = "Arial"
LABEL_SIZE = 18
INSET      = 4   # pt


def pt(x):
    return x * PT_PER_IN


def add_raster(path, x_off, y_off, w_in, h_in):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    panel = Image(pt(w_in), pt(h_in), path)
    panel.move(x_off, y_off)
    return panel


def add_svg(path, x_off, y_off, w_in, native_w_in):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    scale = pt(w_in) / pt(native_w_in)
    panel = SVG(path)
    panel.scale(scale)
    panel.move(x_off, y_off)
    return panel


def stitch(scenario, out_path):
    fig_dir = f"figs/{scenario}"

    total_w = pt(TARGET_W_IN)
    total_h = pt(H_A_IN + H_B_IN)

    row1_y = 0
    row2_y = pt(H_A_IN)

    elements = [Rect(total_w, total_h, fill="white")]

    elements.append(add_raster(
        f"{fig_dir}/RR_maps/patched_maps.jpg",
        0, row1_y, TARGET_W_IN, H_A_IN
    ))
    elements.append(add_svg(
        f"{fig_dir}/state_heatmap_clustered.svg",
        0, row2_y, TARGET_W_IN, NATIVE_W_B_IN
    ))

    for letter, lx, ly in [
        ("A", INSET, row1_y + INSET),
        ("B", INSET, row2_y + INSET),
    ]:
        elements.append(
            Text(letter, lx, ly + LABEL_SIZE,
                 size=LABEL_SIZE, weight="bold", font=LABEL_FONT)
        )

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Figure(f"{total_w}pt", f"{total_h}pt", *elements).save(out_path)
    total_h_in = H_A_IN + H_B_IN
    print(f"Saved: {out_path}  ({TARGET_W_IN:.3f} x {total_h_in:.3f} in)")
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
    parser.add_argument("--out",  default="manuscript/figures/geo_clust.svg")
    parser.add_argument("--pdf",  action="store_true")
    parser.add_argument("--pdf-out", default=None)
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
