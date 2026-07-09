#!/usr/bin/env python3
"""Build a submission DOCX with manuscript text, legends, and embedded figures.

This is a formatting/export helper. It reads the curated manuscript and figure
legend files, appends each figure image after its legend, then calls pandoc to
embed the images in a DOCX file.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path


ROOT = Path(os.environ.get("PROJ_ROOT", ".")).resolve()
MANUSCRIPT = ROOT / "docs" / "VS_manuscript_submission.md"
LEGENDS = ROOT / "results" / "figures_pub" / "figure_legends.md"
OUT_MD = ROOT / "docs" / "VS_manuscript_submission_BMC_embedded.md"
OUT_DOCX = ROOT / "docs" / "VS_manuscript_submission_BMC_embedded.docx"


FIGURE_IMAGES = {
    "Figure 1": ROOT / "results" / "figures_submission" / "Fig1.png",
    "Figure 2": ROOT / "results" / "figures_submission" / "Fig2.png",
    "Figure 3": ROOT / "results" / "figures_submission" / "Fig3.png",
    "Figure 4": ROOT / "results" / "figures_submission" / "Fig4.png",
    "Figure 5": ROOT / "results" / "figures_submission" / "Fig5.png",
    "Supplementary Figure S1": ROOT / "results" / "figures_submission" / "FigS1.png",
    "Supplementary Figure S2": ROOT / "results" / "figures_submission" / "FigS2.png",
    "Supplementary Figure S3": ROOT / "results" / "figures_submission" / "FigS3.png",
    "Supplementary Figure S4": ROOT / "results" / "figures_submission" / "FigS4.png",
    "Supplementary Figure S5": ROOT / "results" / "figures_submission" / "FigS5.png",
    "Supplementary Figure S6": ROOT / "results" / "figures_submission" / "FigS6.png",
}


def strip_existing_figure_note(text: str) -> str:
    marker = "\n## Figure legends\n"
    if marker not in text:
        return text.rstrip() + "\n"
    return text.split(marker, 1)[0].rstrip() + "\n"


def legend_text_with_images(text: str) -> str:
    """Insert each figure image immediately after its legend block."""
    out_lines: list[str] = []
    current_title: str | None = None
    inserted: set[str] = set()

    def flush_current_image() -> None:
        nonlocal current_title
        if current_title is None:
            return
        image = FIGURE_IMAGES[current_title]
        out_lines.append("")
        image_ref = image.relative_to(ROOT).as_posix()
        out_lines.append(f"![{current_title}]({image_ref}){{width=6.5in}}")
        out_lines.append("")
        inserted.add(current_title)
        current_title = None

    for line in text.splitlines():
        figure_match = re.match(r"^### (Figure \d+|Supplementary Figure S\d+)\.", line)
        heading_match = re.match(r"^#{2,3} ", line)
        if heading_match and current_title is not None:
            flush_current_image()
        out_lines.append(line)
        if figure_match:
            current_title = figure_match.group(1)
    flush_current_image()

    missing = set(FIGURE_IMAGES) - inserted
    extra = inserted - set(FIGURE_IMAGES)
    if missing or extra:
        details = []
        if missing:
            details.append("Missing legend blocks: " + ", ".join(sorted(missing)))
        if extra:
            details.append("Unexpected legend blocks: " + ", ".join(sorted(extra)))
        raise SystemExit("\n".join(details))

    return "\n".join(out_lines).rstrip() + "\n"


def apply_bmc_docx_formatting(docx_path: Path) -> None:
    """Apply BMC-friendly DOCX details that pandoc does not emit by default."""
    footer_rel_id = "rIdBMCFooter1"
    footer_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <w:p>
    <w:pPr><w:jc w:val="center" /></w:pPr>
    <w:r><w:fldChar w:fldCharType="begin" /></w:r>
    <w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>
    <w:r><w:fldChar w:fldCharType="separate" /></w:r>
    <w:r><w:t>1</w:t></w:r>
    <w:r><w:fldChar w:fldCharType="end" /></w:r>
  </w:p>
</w:ftr>
"""

    with zipfile.ZipFile(docx_path, "r") as zin:
        entries = {name: zin.read(name) for name in zin.namelist()}

    document_xml = entries["word/document.xml"].decode("utf-8")
    document_xml = re.sub(
        r'(<pic:cNvPr descr=")(?:[^"]*/)?(FigS?\d+\.png")',
        r"\1\2",
        document_xml,
    )
    if footer_rel_id not in document_xml:
        insert = (
            f'<w:footerReference w:type="default" r:id="{footer_rel_id}" />'
            '<w:lnNumType w:countBy="1" />'
        )
        document_xml = re.sub(r"(<w:sectPr\b[^>]*>)", r"\1" + insert, document_xml, count=1)
    elif "w:lnNumType" not in document_xml:
        document_xml = re.sub(
            r"(<w:sectPr\b[^>]*>)",
            r'\1<w:lnNumType w:countBy="1" />',
            document_xml,
            count=1,
        )
    entries["word/document.xml"] = document_xml.encode("utf-8")

    styles_xml = entries["word/styles.xml"].decode("utf-8")
    styles_xml = styles_xml.replace(
        '<w:spacing w:after="180" w:before="180" />',
        '<w:spacing w:after="180" w:before="180" w:line="480" w:lineRule="auto" />',
    )
    entries["word/styles.xml"] = styles_xml.encode("utf-8")

    rels_xml = entries["word/_rels/document.xml.rels"].decode("utf-8")
    if footer_rel_id not in rels_xml:
        rels_xml = rels_xml.replace(
            "</Relationships>",
            f'<Relationship Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
            f'Id="{footer_rel_id}" Target="footer1.xml" /></Relationships>',
        )
    entries["word/_rels/document.xml.rels"] = rels_xml.encode("utf-8")

    content_types_xml = entries["[Content_Types].xml"].decode("utf-8")
    if "/word/footer1.xml" not in content_types_xml:
        content_types_xml = content_types_xml.replace(
            "</Types>",
            '<Override PartName="/word/footer1.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml" /></Types>',
        )
    entries["[Content_Types].xml"] = content_types_xml.encode("utf-8")
    entries["word/footer1.xml"] = footer_xml.encode("utf-8")

    with tempfile.NamedTemporaryFile(delete=False, dir=docx_path.parent, suffix=".docx") as tmp:
        tmp_path = Path(tmp.name)
    try:
        with zipfile.ZipFile(tmp_path, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for name, content in entries.items():
                zout.writestr(name, content)
        shutil.move(str(tmp_path), str(docx_path))
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def main() -> None:
    missing = [str(path) for path in FIGURE_IMAGES.values() if not path.exists()]
    if missing:
        raise SystemExit("Missing figure images:\n" + "\n".join(missing))

    manuscript_text = strip_existing_figure_note(MANUSCRIPT.read_text())
    legends_text = legend_text_with_images(LEGENDS.read_text())

    OUT_MD.write_text(manuscript_text + "\n## Figure legends and embedded figures\n\n" + legends_text)
    subprocess.run(
        ["pandoc", str(OUT_MD), "-o", str(OUT_DOCX), "--resource-path", str(ROOT)],
        check=True,
        cwd=str(ROOT),
    )
    apply_bmc_docx_formatting(OUT_DOCX)
    print(f"Wrote {OUT_MD}")
    print(f"Wrote {OUT_DOCX}")


if __name__ == "__main__":
    main()
