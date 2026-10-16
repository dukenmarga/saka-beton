from __future__ import annotations  # 1. Must be the first line
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biaxial_analysis import ConcreteInteractionDiagram

from dataclasses import dataclass
import os

from fpdf import FPDF
from fpdf.fonts import TextStyle
from markdown_it import MarkdownIt

import geometry



@dataclass
class Report:
    def __init__(self, parent: ConcreteInteractionDiagram) -> None:
        self.parent = parent
        self.images: dict[str, str] = {}

    def prepare_report(self):
        self.markdown_text = f"""
# Biaxial Interaction Diagram Report

## Material

{self._table(self._material_rows())}
"""

        self.markdown_text += f"""
## Dimension

{self._table(self._dimension_rows())}

{self._figure("section", "Concrete Section", 170)}
"""

        self.markdown_text += f"""
## Reinforcement

{self._table(self._reinforcement_rows())}
"""

        self.markdown_text += f"""
## Analysis

{self._loads_table()}

{self._loads_verdict()}

{self._figure("diagram", "PM Interaction Diagram and Ultimate Loads", 470)}
"""

    def _table(self, rows: list[tuple[str, str]]) -> str:
        """
        Render label/value pairs as a borderless two column table. Column widths
        belong on the first row only: fpdf2 ignores them elsewhere and warns.
        """
        cells = []
        for i, (label, value) in enumerate(rows):
            w = (' width="30%"', ' width="2%"', ' width="68%"') if i == 0 else ("", "", "")
            cells.append(
                f"""<tr>
    <td{w[0]}><strong>{label}</strong></td>
    <td{w[1]}>:</td>
    <td{w[2]}>{value}</td>
</tr>"""
            )
        joined = "\n".join(cells)
        return f'<table border="0" width="100%" padding="0px">\n{joined}\n</table>'

    def _material_rows(self) -> list[tuple[str, str]]:
        material = self.parent.material
        return [
            ("fc'", f"{material.fc_:g} MPa"),
            ("fy", f"{material.fy:g} MPa"),
            ("Transversal rebar", f"{self.parent.sni.transverse}"),
        ]

    def _dimension_rows(self) -> list[tuple[str, str]]:
        section = self.parent.section
        shape = section.get("shape", "")

        if shape == "rectangle":
            rows = [
                ("Shape", "Rectangle"),
                ("Width", f"{section['width']:g} mm"),
                ("Height", f"{section['height']:g} mm"),
            ]
        elif shape == "circle":
            rows = [
                ("Shape", "Circle"),
                ("Diameter", f"{section['diameter']:g} mm"),
            ]
        elif shape == "hollow-rectangle":
            rows = [
                ("Shape", "Hollow rectangle"),
                ("Width", f"{section['width']:g} mm"),
                ("Height", f"{section['height']:g} mm"),
                ("Wall thickness", f"{section['thickness']:g} mm"),
            ]
        elif shape == "hollow-circle":
            rows = [
                ("Shape", "Hollow circle"),
                ("Outer diameter", f"{section['outer_diameter']:g} mm"),
                ("Inner diameter", f"{section['inner_diameter']:g} mm"),
            ]
        elif shape.startswith("L-shape") or shape.startswith("T-shape"):
            rows = [
                ("Shape", shape),
                ("A", f"{section['a']:g} mm"),
                ("B", f"{section['b']:g} mm"),
                ("C", f"{section['c']:g} mm"),
                ("D", f"{section['d']:g} mm"),
            ]
        else:
            # No shape was chosen from the form: the section was built straight
            # from coordinates, so only the boundary itself can be described.
            x = self.parent.concrete_boundary["x"]
            y = self.parent.concrete_boundary["y"]
            rows = [
                ("Shape", "Arbitrary (from coordinates)"),
                ("Bounding width", f"{x.max() - x.min():g} mm"),
                ("Bounding height", f"{y.max() - y.min():g} mm"),
            ]

        rows.append(("Overall area", f"{self.parent.outer_concrete_area:,.0f} mm²"))
        rows.append(("Opening area", f"{self.parent.opening_area:,.0f} mm²"))
        rows.append(("Gross area (Ag)", f"{self.parent.gross_concrete_area:,.0f} mm²"))
        return rows

    def _reinforcement_rows(self) -> list[tuple[str, str]]:
        section = self.parent.section
        rebar = self.parent.rebar
        diameters = rebar["diameter"].value_counts().sort_index()
        rebar_area = float(
            geometry.calculate_rebar_area(rebar["diameter"].to_numpy(copy=True)).sum()
        )
        gross_area = self.parent.gross_concrete_area

        rows = [
            ("Number of rebar", f"{len(rebar)}"),
            (
                "Rebar",
                ", ".join(f"{n} \u00d7 D{dia:g}" for dia, n in diameters.items()),
            ),
        ]
        if "cover" in section:
            rows.append(("Cover", f"{section['cover']:g} mm"))
        rows.append(("Rebar area", f"{rebar_area:,.0f} mm²"))
        rows.append(("Net concrete area", f"{gross_area - rebar_area:,.0f} mm²"))
        rows.append(("Rebar ratio", f"{100 * rebar_area / gross_area:.2f} %"))
        return rows

    def _figure(self, name: str, heading: str, width: int) -> str:
        """Heading and image for one figure saved by Plot.save_for_report()."""
        if name not in self.images:
            return ""
        return f'### {heading}\n\n<img src="{self.images[name]}" width="{width}">'

    def _loads_table(self) -> str:
        """One row per load combination, with the check against phi capacity."""
        loads = self.parent.loads
        head = """<tr>
    <th width="10%">No</th>
    <th width="23%">Pu (kN)</th>
    <th width="23%">Mux (kN.m)</th>
    <th width="23%">Muy (kN.m)</th>
    <th width="21%">Status</th>
</tr>"""
        rows = "\n".join(
            f"""<tr>
    <td>{i + 1}</td>
    <td>{Pu:g}</td>
    <td>{Mux:g}</td>
    <td>{Muy:g}</td>
    <td>{"Safe" if safe else "Unsafe"}</td>
</tr>"""
            for i, (Pu, Mux, Muy, safe) in enumerate(
                zip(loads["Pu"], loads["Mux"], loads["Muy"], loads["safe"])
            )
        )
        return f'<table border="1" width="100%">\n{head}\n{rows}\n</table>'

    def _loads_verdict(self) -> str:
        """One sentence on whether the section carries every load combination."""
        safe = self.parent.loads["safe"]
        unsafe = len(safe) - sum(safe)
        if unsafe == 0:
            return f"All {len(safe)} load combinations are within the design capacity."
        return (
            f"{unsafe} of {len(safe)} load combinations exceed the design capacity."
        )

    def _use_unicode_font(self, pdf: FPDF) -> None:
        """
        The built in fonts have no glyph for mm², so use the DejaVu font that
        ships with matplotlib. Imported here to keep matplotlib out of the
        workers that import this module.
        """
        import matplotlib

        ttf = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
        pdf.add_font("DejaVu", "", os.path.join(ttf, "DejaVuSans.ttf"))
        pdf.add_font("DejaVu", "B", os.path.join(ttf, "DejaVuSans-Bold.ttf"))
        pdf.set_font("DejaVu", size=11)

    def convert_report(self):
        md = MarkdownIt()
        html_content = md.render(self.markdown_text)

        pdf = FPDF()
        pdf.add_page()
        self._use_unicode_font(pdf)
        pdf.set_text_color(0, 0, 0)

        # Menggunakan fitur internal FPDF2 untuk mencetak HTML
        pdf.write_html(
            html_content,
            font_family="DejaVu",
            tag_styles={
                # Default heading margins are ~5mm on top; trimmed here because
                # fpdf2 already breaks a line after every </table>.
                "h1": TextStyle(font_size_pt=18, t_margin=0, b_margin=1),
                "h2": TextStyle(font_size_pt=14, t_margin=0, b_margin=0.5),
                "h3": TextStyle(font_size_pt=12, t_margin=0, b_margin=0.5),
            },
        )

        # 4. Output menjadi file PDF
        pdf.output("analysis_report.pdf")