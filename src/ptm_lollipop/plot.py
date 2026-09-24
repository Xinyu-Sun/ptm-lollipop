from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path

import fitz
from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

from .models import ProteinRecord

PTM_COLORS = {
    "phosphorylation": "#D62728",
    "sumoylation": "#9467BD",
    "ubiquitination": "#2CA02C",
    "acetylation": "#17BECF",
    "methylation": "#E377C2",
    "succinylation": "#8C564B",
    "glutathionylation": "#7F7F7F",
    "other": "#7F7F7F",
}


def slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower() or "proteins"


def plot_all_categories(records: list[ProteinRecord], outdir: str | Path, dpi: int = 300) -> list[Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    categories = []
    for record in records:
        if record.category not in categories:
            categories.append(record.category)

    written: list[Path] = []
    for category in categories:
        stem = slugify(category)
        pdf_path = out / f"{stem}_ptm_lollipop.pdf"
        png_path = out / f"{stem}_ptm_lollipop.png"
        plot_category([r for r in records if r.category == category], category, pdf_path, png_path, dpi=dpi)
        written.extend([pdf_path, png_path])
    return written


def plot_category(records: list[ProteinRecord], category: str, path_pdf: Path, path_png: Path, dpi: int = 300) -> None:
    if not records:
        return
    recs = sorted(records, key=lambda record: record.gene)
    max_len = max(record.length for record in recs)
    width_pt = letter[0]
    row_h = 23
    legend_step = 10
    legend_height = 18 + max(1, len({ptm.family for record in recs for ptm in record.ptms})) * legend_step
    top = max(108, 70 + legend_height)
    bottom = 52
    height_pt = max(210, top + len(recs) * row_h + bottom)
    left_label = 34
    x0 = 230
    x1 = width_pt - 42
    aa_scale = (x1 - x0) / max_len

    families = sorted(
        {ptm.family for record in recs for ptm in record.ptms},
        key=lambda family: list(PTM_COLORS).index(family) if family in PTM_COLORS else 99,
    )

    def aa_to_x(pos: int) -> float:
        return x0 + (pos - 1) * aa_scale

    def draw_pdf(c) -> None:
        def s(value: float) -> float:
            return value

        def color(hex_color: str):
            return HexColor(hex_color)

        c.setFillColor(color("#202124"))
        c.setFont("Helvetica-Bold", 12)
        c.drawString(s(x0), s(height_pt - 27), category)
        c.setFont("Helvetica", 7)
        c.setFillColor(color("#5f6368"))
        c.drawRightString(s(x1), s(18), f"{max_len:,} aa")
        c.setStrokeColor(color("#9aa0a6"))
        c.setLineWidth(s(0.7))
        c.line(s(x0), s(30), s(x1), s(30))

        tick_max = int(math.ceil(max_len / 500.0) * 500)
        for tick in range(500, tick_max + 1, 500):
            if tick > max_len:
                continue
            x = aa_to_x(tick)
            show_tick_label = x1 - x >= 58
            c.setStrokeColor(color("#9aa0a6"))
            c.line(s(x), s(26), s(x), s(34))
            if show_tick_label:
                c.setFont("Helvetica", 6.5)
                c.setFillColor(color("#5f6368"))
                c.drawCentredString(s(x), s(14), str(tick))

        for i, record in enumerate(recs):
            row_y = height_pt - top - i * row_h
            c.setFillColor(color("#202124"))
            c.setFont("Helvetica-Bold", 8.5)
            c.drawString(s(left_label), s(row_y - 2), record.gene)
            c.setFillColor(color("#5f6368"))
            c.setFont("Helvetica", 6.5)
            c.drawString(s(left_label + 78), s(row_y - 2), f"{record.systematic} | {record.length} aa")
            c.setStrokeColor(color("#343a40"))
            c.setLineWidth(s(1.7))
            c.line(s(x0), s(row_y), s(aa_to_x(record.length)), s(row_y))

            for start, end in record.disorder:
                xa, xb = aa_to_x(start), aa_to_x(end)
                c.setFillColor(color("#5B5B5B"))
                c.rect(s(xa), s(row_y - 3.3), s(max(1.0, xb - xa)), s(6.6), stroke=0, fill=1)

            sites_by_position = defaultdict(list)
            for ptm in record.ptms:
                sites_by_position[ptm.site].append(ptm)
            for site, ptms in sites_by_position.items():
                for j, ptm in enumerate(ptms):
                    ptm_color = color(PTM_COLORS.get(ptm.family, PTM_COLORS["other"]))
                    x = aa_to_x(site)
                    offset = (j - (len(ptms) - 1) / 2) * 3.0
                    c.setStrokeColor(ptm_color)
                    c.setLineWidth(s(0.85))
                    c.line(s(x), s(row_y + 4.5), s(x), s(row_y + 11.5))
                    c.setFillColor(ptm_color)
                    c.circle(s(x), s(row_y + 14 + offset), s(2.0), stroke=0, fill=1)

        legend_x = width_pt - 193
        legend_y = height_pt - 29
        legend_font_size = 7 if len(families) <= 6 else 6.2
        legend_marker = 2.2 if len(families) <= 6 else 1.8
        legend_step_dynamic = 10 if len(families) <= 6 else 8.5
        c.setFont("Helvetica", legend_font_size)
        c.setFillColor(color("#5B5B5B"))
        c.rect(s(legend_x), s(legend_y - 3.5), s(10), s(4.5), stroke=0, fill=1)
        c.setFillColor(color("#202124"))
        c.drawString(s(legend_x + 14), s(legend_y - 3.5), "MobiDB-lite disorder")
        for k, family in enumerate(families):
            ly = legend_y - 12 - k * legend_step_dynamic
            c.setFillColor(color(PTM_COLORS.get(family, PTM_COLORS["other"])))
            c.circle(s(legend_x + 5), s(ly), s(legend_marker), stroke=0, fill=1)
            c.setFillColor(color("#202124"))
            c.drawString(s(legend_x + 14), s(ly - 2.2), family.capitalize())

    pdf = canvas.Canvas(str(path_pdf), pagesize=(width_pt, height_pt))
    draw_pdf(pdf)
    pdf.save()

    render_pdf_to_png(path_pdf, path_png, dpi=dpi)


def render_pdf_to_png(path_pdf: Path, path_png: Path, dpi: int = 300) -> None:
    zoom = dpi / 72
    with fitz.open(path_pdf) as doc:
        page = doc.load_page(0)
        pixmap = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
        pixmap.save(path_png)
