from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
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
    top = 108
    bottom = 52
    height_pt = max(210, top + len(recs) * row_h + bottom)
    left_label = 34
    x0 = 230
    x1 = width_pt - 42
    scale = (x1 - x0) / max_len

    families = sorted(
        {ptm.family for record in recs for ptm in record.ptms},
        key=lambda family: list(PTM_COLORS).index(family) if family in PTM_COLORS else 99,
    )

    def aa_to_x(pos: int) -> float:
        return x0 + (pos - 1) * scale

    def draw_common(drawer, is_pdf: bool, dpi_scale: float = 1.0) -> None:
        def s(value: float) -> float:
            return value * dpi_scale

        def y(value: float) -> float:
            return s(height_pt - value)

        def color(hex_color: str):
            if is_pdf:
                return HexColor(hex_color)
            cleaned = hex_color.lstrip("#")
            return tuple(int(cleaned[i:i + 2], 16) for i in (0, 2, 4))

        if is_pdf:
            c = drawer
            c.setFillColor(color("#202124"))
            c.setFont("Helvetica-Bold", 12)
            c.drawString(s(x0), s(height_pt - 27), category)
            c.setFont("Helvetica", 7)
            c.setFillColor(color("#5f6368"))
            c.drawRightString(s(x1), s(18), f"{max_len:,} aa")
            c.setStrokeColor(color("#9aa0a6"))
            c.setLineWidth(s(0.7))
            c.line(s(x0), s(30), s(x1), s(30))
        else:
            d = drawer
            fonts = _fonts()
            d.rectangle([0, 0, s(width_pt), s(height_pt)], fill="white")
            d.text((s(x0), s(18)), category, font=fonts["title"], fill=color("#202124"))
            d.text((s(x1 - 48), y(18)), f"{max_len:,} aa", font=fonts["small"], fill=color("#5f6368"))
            d.line([s(x0), y(30), s(x1), y(30)], fill=color("#9aa0a6"), width=max(1, int(s(0.7))))

        tick_max = int(math.ceil(max_len / 500.0) * 500)
        for tick in range(500, tick_max + 1, 500):
            if tick > max_len:
                continue
            x = aa_to_x(tick)
            if is_pdf:
                c.setStrokeColor(color("#9aa0a6"))
                c.line(s(x), s(26), s(x), s(34))
                c.setFont("Helvetica", 6.5)
                c.setFillColor(color("#5f6368"))
                c.drawCentredString(s(x), s(14), str(tick))
            else:
                d.line([s(x), y(34), s(x), y(26)], fill=color("#9aa0a6"), width=max(1, int(s(0.7))))
                d.text((s(x - 8), y(14)), str(tick), font=fonts["small"], fill=color("#5f6368"))

        for i, record in enumerate(recs):
            row_y = height_pt - top - i * row_h
            if is_pdf:
                c.setFillColor(color("#202124"))
                c.setFont("Helvetica-Bold", 8.5)
                c.drawString(s(left_label), s(row_y - 2), record.gene)
                c.setFillColor(color("#5f6368"))
                c.setFont("Helvetica", 6.5)
                c.drawString(s(left_label + 78), s(row_y - 2), f"{record.systematic} | {record.length} aa")
                c.setStrokeColor(color("#343a40"))
                c.setLineWidth(s(1.7))
                c.line(s(x0), s(row_y), s(aa_to_x(record.length)), s(row_y))
            else:
                d.text((s(left_label), y(row_y + 8)), record.gene, font=fonts["bold"], fill=color("#202124"))
                d.text((s(left_label + 78), y(row_y + 7)), f"{record.systematic} | {record.length} aa", font=fonts["small"], fill=color("#5f6368"))
                d.line([s(x0), y(row_y), s(aa_to_x(record.length)), y(row_y)], fill=color("#343a40"), width=max(2, int(s(1.7))))

            for start, end in record.disorder:
                xa, xb = aa_to_x(start), aa_to_x(end)
                if is_pdf:
                    c.setFillColor(color("#5B5B5B"))
                    c.rect(s(xa), s(row_y - 3.3), s(max(1.0, xb - xa)), s(6.6), stroke=0, fill=1)
                else:
                    d.rectangle([s(xa), y(row_y + 4), s(max(xa + 1, xb)), y(row_y - 4)], fill=color("#5B5B5B"))

            sites_by_position = defaultdict(list)
            for ptm in record.ptms:
                sites_by_position[ptm.site].append(ptm)
            for site, ptms in sites_by_position.items():
                for j, ptm in enumerate(ptms):
                    ptm_color = color(PTM_COLORS.get(ptm.family, PTM_COLORS["other"]))
                    x = aa_to_x(site)
                    offset = (j - (len(ptms) - 1) / 2) * 3.0
                    if is_pdf:
                        c.setStrokeColor(ptm_color)
                        c.setLineWidth(s(0.85))
                        c.line(s(x), s(row_y + 4.5), s(x), s(row_y + 11.5))
                        c.setFillColor(ptm_color)
                        c.circle(s(x), s(row_y + 14 + offset), s(2.0), stroke=0, fill=1)
                    else:
                        d.line([s(x), y(row_y + 11.5), s(x), y(row_y + 4.5)], fill=ptm_color, width=max(1, int(s(0.85))))
                        rad = s(2.2)
                        cy = y(row_y + 14 + offset)
                        d.ellipse([s(x) - rad, cy - rad, s(x) + rad, cy + rad], fill=ptm_color, outline="white")

        legend_x = width_pt - 193
        legend_y = height_pt - 29
        if is_pdf:
            c.setFont("Helvetica", 7)
            c.setFillColor(color("#5B5B5B"))
            c.rect(s(legend_x), s(legend_y - 3.5), s(10), s(4.5), stroke=0, fill=1)
            c.setFillColor(color("#202124"))
            c.drawString(s(legend_x + 14), s(legend_y - 3.5), "MobiDB-lite disorder")
            for k, family in enumerate(families):
                ly = legend_y - 12 - k * 10
                c.setFillColor(color(PTM_COLORS.get(family, PTM_COLORS["other"])))
                c.circle(s(legend_x + 5), s(ly), s(2.2), stroke=0, fill=1)
                c.setFillColor(color("#202124"))
                c.drawString(s(legend_x + 14), s(ly - 2.2), family.capitalize())
        else:
            d.rectangle([s(legend_x), s(24), s(legend_x + 10), s(28.5)], fill=color("#5B5B5B"))
            d.text((s(legend_x + 14), s(19)), "MobiDB-lite disorder", font=fonts["small"], fill=color("#202124"))
            for k, family in enumerate(families):
                ly = 36 + k * 10
                ptm_color = color(PTM_COLORS.get(family, PTM_COLORS["other"]))
                d.ellipse([s(legend_x + 3.3), s(ly - 2.2), s(legend_x + 6.7), s(ly + 2.2)], fill=ptm_color)
                d.text((s(legend_x + 14), s(ly - 4.4)), family.capitalize(), font=fonts["small"], fill=color("#202124"))

    pdf = canvas.Canvas(str(path_pdf), pagesize=(width_pt, height_pt))
    draw_common(pdf, is_pdf=True)
    pdf.save()

    scale = dpi / 72
    image = Image.new("RGB", (int(width_pt * scale), int(height_pt * scale)), "white")
    draw = ImageDraw.Draw(image)
    draw_common(draw, is_pdf=False, dpi_scale=scale)
    image.save(path_png, dpi=(dpi, dpi))


def _fonts() -> dict[str, ImageFont.FreeTypeFont | ImageFont.ImageFont]:
    return {
        "small": _font(28),
        "bold": _font(36, bold=True),
        "title": _font(44, bold=True),
    }


def _font(size: int, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for path in candidates:
        if path and Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()
