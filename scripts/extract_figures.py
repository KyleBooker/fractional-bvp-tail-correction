"""One-shot helper to extract paper figures from GMp.pdf into docs/figures/.

Not part of the published code — just a utility we used once to render
README assets from the paper PDF. Safe to delete or keep for reproducibility.
"""
import fitz
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "GMp.pdf"
OUT = ROOT / "docs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

DPI = 220
ZOOM = DPI / 72.0
MAT = fitz.Matrix(ZOOM, ZOOM)

# Trim page-edge whitespace and headers/footers.
LEFT_PT, RIGHT_PT = 78, 522

# (page_index_0based, top_pt, bottom_pt, output_name, caption_for_log)
TARGETS = [
    (3,  120, 440, "fig01_tail_correction_effect.png",
     "Fig 1: spike with vs. without tail correction (the headline figure)"),
    (8,  115, 495, "fig06_relative_error_schemes.png",
     "Fig 6: relative error across the four tail-condition schemes"),
    (10, 115, 740, "fig07_solution_shapes_p2.png",
     "Fig 7: solution u(y;2,gamma) for several gamma values, core + tail"),
    (14, 115, 490, "fig11_spike_peak_vs_gamma.png",
     "Fig 11: spike peak u(0;p,gamma) vs gamma for p in [2,5]"),
]


def main():
    doc = fitz.open(PDF)
    for page_idx, top, bottom, name, caption in TARGETS:
        page = doc[page_idx]
        clip = fitz.Rect(LEFT_PT, top, RIGHT_PT, bottom)
        pix = page.get_pixmap(matrix=MAT, clip=clip, alpha=False)
        out_path = OUT / name
        pix.save(str(out_path))
        print(f"  {out_path.name:42s}  {caption}")
    doc.close()


if __name__ == "__main__":
    main()
