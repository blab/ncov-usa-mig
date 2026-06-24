#!/usr/bin/env python3
"""Stitch the age-analysis figures into a single multi-panel SVG.

Layout (subfigure letters A-D; contact correlation moved to a supplement):
  +------------------------------------------------------+
  |             A | full heatmap (expanded, centered)    |
  +------------------------------------------------------+
  | B | young_0_24 heatmap | C | subsets_row             |
  |                        | D | deviance                |
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

# ----- Panel A (top row): full heatmap, expanded and centered -----
# full.svg is saved on a 6 (W) x 5.3 (H) canvas (trimmed to remove vertical
# whitespace around the heatmap; keep A_SRC_H_IN in sync with the ggsave height
# in age_heatmap.R). A_SIDE_IN is the single knob for how wide A renders; its
# height follows from the native aspect ratio.
A_SRC_W_IN = 6
A_SRC_H_IN = 5.3
A_SIDE_IN = 10                              # rendered width of panel A
A_SCALE = A_SIDE_IN / A_SRC_W_IN
A_H_IN = A_SIDE_IN * (A_SRC_H_IN / A_SRC_W_IN)
A_X_IN = (PANEL_W_IN - A_SIDE_IN) / 2.0     # horizontal centering offset

# ----- Bottom row: young heatmap (B) left, subsets(C)/deviance(D) stack right -----
# Mirrors the old top-right composite: subsets row sits on top of the deviance
# line plot. The young heatmap is a square filling the row height on the left.
RIGHT_W_IN = 5             # width of the subsets/deviance stack (native subsets width)
BOTTOM_H_IN = 5            # subsets (2") + deviance (3")
SUBSETS_H_IN = 2.0

YOUNG_SRC_IN = 6           # young_0_24.svg native canvas (square)
YOUNG_SIDE_IN = BOTTOM_H_IN
YOUNG_SCALE = YOUNG_SIDE_IN / YOUNG_SRC_IN  # square heatmap fills the row height

# Deviance saved at 5 x 3.33; scaled to its reduced height allotment and nudged
# right so its content clears the "D" letter label.
DEV_NATIVE_H_IN = 3.33
DEV_INSET_IN = 0.25
DEV_H_IN = BOTTOM_H_IN - SUBSETS_H_IN  # 3.0
DEV_SCALE = DEV_H_IN / DEV_NATIVE_H_IN  # height-driven shrink

TOTAL_H_IN = A_H_IN + BOTTOM_H_IN


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

    # ----- Panel A (top row): full heatmap, expanded and centered -----
    a_x = LEFT_MARGIN + in_pt(A_X_IN)
    elements.append(
        add_panel(f"{fig_dir}/age_heatmaps/full.svg",
                  a_x, 0, scale=A_SCALE)
    )

    # ----- Bottom row -----
    bc_y = in_pt(A_H_IN)
    right_x = LEFT_MARGIN + in_pt(YOUNG_SIDE_IN)

    # Panel B: young 0-24yo heatmap (bottom-left square).
    elements.append(
        add_panel(f"{fig_dir}/age_heatmaps/young_0_24.svg",
                  LEFT_MARGIN, bc_y, scale=YOUNG_SCALE)
    )
    # Panel C: subsets row (top of the bottom-right stack), native 5 x 2.
    elements.append(
        add_panel(f"{fig_dir}/age_heatmaps/subsets_row.svg",
                  right_x, bc_y)
    )
    # Panel D: deviance (below subsets), nudged right so it clears the "D" label.
    elements.append(
        add_panel(f"{fig_dir}/age_RR_deviance_geographic_compact.svg",
                  right_x + in_pt(DEV_INSET_IN),
                  bc_y + in_pt(SUBSETS_H_IN),
                  scale=DEV_SCALE)
    )

    # ----- Subfigure letters -----
    # A: full heatmap (top-left edge of the centered panel)
    # B: young 0-24 heatmap (left gutter, bottom row)
    # C: subsets_row top-left (inside, top of right stack)
    # D: deviance top-left (inside, below subsets row)
    INSET = 4
    letters = [
        ("A", a_x,                       INSET),
        ("B", LETTER_X,                  bc_y),
        ("C", right_x + INSET,           bc_y + INSET),
        ("D", right_x + INSET,           bc_y + in_pt(SUBSETS_H_IN) + INSET),
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
