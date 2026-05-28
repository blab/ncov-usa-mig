#!/usr/bin/env python3
"""Stitch the age-analysis figures into a single multi-panel SVG.

Layout (subfigure letters A-E):
  +------------------------------------------------------+
  | A | full heatmap         | B | subsets_row          |
  |   |                      |   | C | deviance         |
  +------------------------------------------------------+
  | D | contact correlation  | E | vaccine nRR          |
  +------------------------------------------------------+
"""

import argparse
import os

from svgutils.compose import Figure, SVG, Image, Text
from svg_helpers import Rect



PT_PER_IN = 72
PANEL_W_IN = 10            # full width of the figure (sans gutter)
W = PANEL_W_IN * PT_PER_IN

LABEL_FONT = "Arial"
LABEL_SIZE = 18
LEFT_MARGIN = 24
LETTER_X = 8

# Row heights (inches)
A_H_IN = 5
BC_H_IN = 3
TOTAL_H_IN = A_H_IN + BC_H_IN

# Panel A internal layout (within 10" wide, 5" tall)
A_FULL_W_IN = 5            # left half: full heatmap area
A_RIGHT_W_IN = 5
A_SUBSETS_H_IN = 1.67      # top of right column
A_DEV_H_IN = A_H_IN - A_SUBSETS_H_IN  # 3.33

# Panel B/C share the bottom row (5" wide each)
BC_PANEL_W_IN = 5


def in_pt(x):
    return x * PT_PER_IN


def add_panel(path, x_off, y_off, scale=None, w_in=None, h_in=None):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    ext = os.path.splitext(path)[1].lower()
    if ext == ".svg":
        panel = SVG(path)
        if scale is not None:
            panel.scale(scale)
        panel.move(x_off, y_off)
        return panel

    if ext in (".jpg", ".jpeg", ".png"):
        if w_in is None or h_in is None:
            raise ValueError(f"raster {path} needs w_in and h_in")
        panel = Image(w_in * PT_PER_IN, h_in * PT_PER_IN, path)
        panel.move(x_off, y_off)
        return panel

    raise ValueError(f"Unsupported file type: {path}")


def stitch(scenario, out_path):
    fig_dir = f"figs/{scenario}"

    total_w = W + LEFT_MARGIN
    total_h = in_pt(TOTAL_H_IN)

    elements = [Rect(total_w, total_h, fill="white")]

    # ----- Panel A: heatmap composite -----
    a_y = 0
    # Full heatmap saved at 6x6 — uniformly scale to fit A_FULL_W_IN square.
    full_scale = A_FULL_W_IN / 6.0
    elements.append(
        add_panel(f"{fig_dir}/age_heatmaps/full.svg",
                  LEFT_MARGIN, a_y, scale=full_scale)
    )
    # Subsets row saved at 5 x 1.67 — slot top-right.
    elements.append(
        add_panel(f"{fig_dir}/age_heatmaps/subsets_row.svg",
                  LEFT_MARGIN + in_pt(A_FULL_W_IN), a_y)
    )
    # Compact deviance saved at 5 x 3.33 — slot below the subset row.
    elements.append(
        add_panel(f"{fig_dir}/age_RR_deviance_geographic_compact.svg",
                  LEFT_MARGIN + in_pt(A_FULL_W_IN),
                  a_y + in_pt(A_SUBSETS_H_IN))
    )

    # ----- Panel B: contact correlation (bottom-left) -----
    bc_y = in_pt(A_H_IN)
    elements.append(
        add_panel(f"{fig_dir}/age_contact_correlation_compact.svg",
                  LEFT_MARGIN, bc_y)
    )

    # ----- Panel C: vaccine period nRR (bottom-right) -----
    elements.append(
        add_panel(f"{fig_dir}/age_time/all_countries/vaccine_period_nRR_fixed_elderly_ci_compact.svg",
                  LEFT_MARGIN + in_pt(BC_PANEL_W_IN), bc_y)
    )

    # ----- Subfigure letters -----
    # A: full heatmap (left gutter, top of panel A)
    # B: subsets_row top-left (inside)
    # C: deviance top-left (inside, below subsets row)
    # D: contact correlation (left gutter, bottom row)
    # E: vaccine top-left (inside, bottom-right panel)
    INSET = 4
    letters = [
        ("A", LETTER_X,                                       0),
        ("B", LEFT_MARGIN + in_pt(A_FULL_W_IN) + INSET,       INSET),
        ("C", LEFT_MARGIN + in_pt(A_FULL_W_IN) + INSET,       in_pt(A_SUBSETS_H_IN) + INSET),
        ("D", LETTER_X,                                       bc_y),
        ("E", LEFT_MARGIN + in_pt(BC_PANEL_W_IN) + INSET,     bc_y + INSET),
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
    parser.add_argument("--out", default="manuscript/figures/age.svg")
    parser.add_argument("--pdf", action="store_true", help="Also export a PDF")
    parser.add_argument("--pdf-out", default=None, help="PDF output path (default: same as --out with .pdf extension)")
    args = parser.parse_args()
    svg_path = stitch(args.scenario, args.out)
    if args.pdf:
        export_pdf(svg_path, args.pdf_out)
