#!/usr/bin/env python3
"""Stitch the region-map + PCoA panels into their own multi-panel SVG.

Layout:
  +----------+-------+
  |    A     |   B   |   bea_region_map.png (raster) / pcoa_V1V2.svg (vector)
  +----------+-------+

(Pulled out of geo_clust.pdf, where they were C/D, since that figure was
getting too tall.)

A's transparent top/bottom padding (an artifact of geom_sf's fixed-aspect
panel expansion) is trimmed before embedding. B (pcoa_V1V2.svg) is exported
with a fixed panel size via save_fixed_panel() in state_time_clust.R, so it
arrives here with no whitespace padding to work around.

Column widths are solved so A and B, scaled to a shared row height, sum to
the full figure width (W_A + W_B = TARGET_W_IN):
    W_A = H / a_aspect,  W_B = H / b_aspect,  W_A + W_B = TARGET_W_IN
    => H * (1/a_aspect + 1/b_aspect) = TARGET_W_IN
"""

import argparse
import os
import re
import subprocess

from PIL import Image as PILImage

from svgutils.compose import Figure, Image, SVG, Text
from svg_helpers import Rect

PT_PER_IN = 72

# Total figure width
TARGET_W_IN = 8.0

LABEL_FONT = "Arial"
LABEL_SIZE = 18
INSET      = 4   # pt


def pt(x):
    return x * PT_PER_IN


def trim_png(path):
    """Crop transparent/whitespace padding from a PNG, caching the result
    alongside the original. Returns (trimmed_path, width_px, height_px)."""
    trimmed_path = os.path.splitext(path)[0] + "_trimmed.png"
    im = PILImage.open(path).convert("RGBA")
    bbox = im.getbbox()
    cropped = im.crop(bbox) if bbox else im
    cropped.save(trimmed_path)
    return trimmed_path, cropped.width, cropped.height


def get_svg_size_in(path):
    """Parse an SVG's width/height (assumed to be in pt) and return inches."""
    with open(path) as f:
        header = f.read(2000)
    w_match = re.search(r"width=['\"]([\d.]+)pt['\"]", header)
    h_match = re.search(r"height=['\"]([\d.]+)pt['\"]", header)
    if not (w_match and h_match):
        raise ValueError(f"Could not parse width/height (pt) from {path}")
    return float(w_match.group(1)) / PT_PER_IN, float(h_match.group(1)) / PT_PER_IN


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

    # Trim A's transparent top/bottom padding, then derive its true aspect ratio
    region_map_path = f"{fig_dir}/bea_region_map.png"
    a_trimmed_path, a_px_w, a_px_h = trim_png(region_map_path)
    a_aspect = a_px_h / a_px_w   # height / width

    # B's native aspect ratio, read directly from the svg (fixed-panel export,
    # no whitespace padding -- see save_fixed_panel() in state_time_clust.R)
    pcoa_path = f"{fig_dir}/clust/pcoa_V1V2.svg"
    b_native_w_in, b_native_h_in = get_svg_size_in(pcoa_path)
    b_aspect = b_native_h_in / b_native_w_in

    # Shared row height, solved so both panels are undistorted and their
    # widths sum to the full figure width
    h_row_in = TARGET_W_IN / (1 / a_aspect + 1 / b_aspect)
    w_a_in = h_row_in / a_aspect
    w_b_in = h_row_in / b_aspect

    total_w = pt(TARGET_W_IN)
    total_h = pt(h_row_in)

    b_x = pt(w_a_in)

    elements = [Rect(total_w, total_h, fill="white")]

    elements.append(add_raster(
        a_trimmed_path,
        0, 0, w_a_in, h_row_in
    ))
    elements.append(add_svg(
        pcoa_path,
        b_x, 0, w_b_in, b_native_w_in
    ))

    for letter, lx, ly in [
        ("A", INSET,      INSET),
        ("B", b_x + INSET, INSET),
    ]:
        elements.append(
            Text(letter, lx, ly + LABEL_SIZE,
                 size=LABEL_SIZE, weight="bold", font=LABEL_FONT)
        )

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Figure(f"{total_w}pt", f"{total_h}pt", *elements).save(out_path)
    print(f"Saved: {out_path}  ({TARGET_W_IN:.3f} x {h_row_in:.3f} in)")
    print(f"  A = {w_a_in:.3f} x {h_row_in:.3f} in, B = {w_b_in:.3f} x {h_row_in:.3f} in")
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
    parser.add_argument("--out",  default="manuscript/figures/geo_pc.svg")
    parser.add_argument("--pdf",  action="store_true")
    parser.add_argument("--pdf-out", default=None)
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
