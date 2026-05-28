#!/usr/bin/env python3
"""Stitch geographic clustering figure panels into a single multi-panel SVG.

Layout:
  +------------------+
  |        A         |   8.000 x 2.133 in  patched_maps.jpg     (raster)
  +----------+-------+
  |    B     |       |   4.290 x 3.754 in  bea_region_map.png   (raster)
  +----------+   D   |
  |    C     |       |   4.290 x 3.575 in  state_heatmap_clustered.svg (vector)
  +----------+-------+
                         3.710 x 7.329 in  pcoa_combined.svg    (vector)

H_B + H_C = H_D = 7.329 in
W_left + W_D = W_A = 8.000 in
"""

import argparse
import os
import subprocess

from svgutils.compose import Figure, Image, SVG, Text
from svg_helpers import Rect

PT_PER_IN = 72

# Total figure width
TARGET_W_IN = 8.0

# Column widths (derived from H_B + H_C = H_D constraint at TARGET_W_IN = 8in)
W_LEFT_IN = 4.290
W_D_IN    = 3.710

# Panel heights
H_A_IN = 2.133   # patched_maps scaled to full width
H_B_IN = 3.754   # region map scaled to W_LEFT_IN
H_C_IN = 3.575   # heatmap scaled to W_LEFT_IN
H_D_IN = 7.329   # H_B + H_C

# Native ggsave widths (used to compute scale factors)
NATIVE_W_A_IN = 15   # patched_maps.jpg
NATIVE_W_B_IN =  8   # bea_region_map.png
NATIVE_W_C_IN = 12   # state_heatmap_clustered.svg
NATIVE_W_D_IN =  3   # pcoa_combined.svg

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
    total_h = pt(H_A_IN + H_D_IN)

    row1_y = 0
    row2_y = pt(H_A_IN)
    row3_y = pt(H_A_IN + H_B_IN)
    d_x    = pt(W_LEFT_IN)

    elements = [Rect(total_w, total_h, fill="white")]

    elements.append(add_raster(
        f"{fig_dir}/RR_maps/patched_maps.jpg",
        0, row1_y, TARGET_W_IN, H_A_IN
    ))
    elements.append(add_raster(
        f"{fig_dir}/bea_region_map.png",
        0, row2_y, W_LEFT_IN, H_B_IN
    ))
    elements.append(add_svg(
        f"{fig_dir}/state_heatmap_clustered.svg",
        0, row3_y, W_LEFT_IN, NATIVE_W_C_IN
    ))
    elements.append(add_svg(
        f"{fig_dir}/clust/pcoa_combined.svg",
        d_x, row2_y, W_D_IN, NATIVE_W_D_IN
    ))

    for letter, lx, ly in [
        ("A", INSET,        row1_y + INSET),
        ("B", INSET,        row2_y + INSET),
        ("C", INSET,        row3_y + INSET),
        ("D", d_x + INSET,  row2_y + INSET),
    ]:
        elements.append(
            Text(letter, lx, ly + LABEL_SIZE,
                 size=LABEL_SIZE, weight="bold", font=LABEL_FONT)
        )

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Figure(f"{total_w}pt", f"{total_h}pt", *elements).save(out_path)
    print(f"Saved: {out_path}  ({TARGET_W_IN:.3f} x {H_A_IN + H_D_IN:.3f} in)")
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
