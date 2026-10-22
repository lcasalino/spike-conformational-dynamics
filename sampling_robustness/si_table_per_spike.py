#!/usr/bin/env python3
# ---------------------------------------------------------------------------
#  Supplementary table: per-spike summary statistics for the 29 RAV spikes.
#
#  For each spike, the value at the first production frame and the mean ± SD
#  over the analysis window (final 352 ns) of the hip, knee and ankle tilt
#  angles and of the NTD triangle area, together with the conformational state
#  of the chain-A RBD. The final row gives the ensemble mean ± SD across the
#  29 spikes. Comparing the t = 0 and window columns shows that each spike's
#  sampled conformation is not predicted by the structure it started from.
#
#  Writes the table in five formats: CSV (data), Markdown, LaTeX (booktabs),
#  HTML and .docx. The .docx is a minimal, dependency-free OOXML package and
#  opens in Word (desktop or online) as a native, editable landscape table.
#
#  Usage:  python3 si_table_per_spike.py
# ---------------------------------------------------------------------------
import os
import csv
import zipfile
import numpy as np

from sampling_robustness import OBS, N_LAST, SPIKES, load_rav, HERE

COLS = ("hip", "knee", "ankle", "ntd")
HEAD = ("Hip, t = 0 (°)", "Hip (°)", "Knee, t = 0 (°)", "Knee (°)",
        "Ankle, t = 0 (°)", "Ankle (°)", "NTD, t = 0 (Å²)", "NTD area (Å²)", "RBD-A")
RBD_UP_CUTOFF = 60.0              # Å; the chain-A RBD-CH distance is bimodal (~53 down, ~68 up)


W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def _esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _cell(text, bold=False, rule=None):
    """One table cell; rule = 'top' | 'bottom' | 'both' draws horizontal rules."""
    b = '<w:b/>' if bold else ''
    edges = {"top": ("top",), "bottom": ("bottom",), "both": ("top", "bottom")}.get(rule, ())
    bd = "".join(f'<w:{e} w:val="single" w:sz="8" w:color="000000"/>' for e in edges)
    bd = f"<w:tcBorders>{bd}</w:tcBorders>" if bd else ""
    return (f'<w:tc><w:tcPr>{bd}</w:tcPr><w:p><w:pPr><w:spacing w:before="20" w:after="20"/>'
            f'<w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:sz w:val="18"/>{b}</w:rPr></w:pPr>'
            f'<w:r><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:sz w:val="18"/>{b}</w:rPr>'
            f'<w:t xml:space="preserve">{_esc(text)}</w:t></w:r></w:p></w:tc>')


def write_docx(path, header, rows, ens_row):
    tr = ["<w:tr>" + "".join(_cell(h, bold=True, rule="both") for h in header) + "</w:tr>"]
    for r in rows:
        tr.append("<w:tr>" + "".join(_cell(c) for c in r) + "</w:tr>")
    tr.append("<w:tr>" + "".join(_cell(c, bold=True, rule="both") for c in ens_row) + "</w:tr>")
    tbl = ('<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/>'
           '<w:tblCellMar><w:left w:w="80" w:type="dxa"/><w:right w:w="80" w:type="dxa"/></w:tblCellMar>'
           '</w:tblPr>' + "".join(tr) + "</w:tbl>")
    doc = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
           f'<w:document xmlns:w="{W}"><w:body>{tbl}'
           f'<w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/>'
           f'<w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134"/></w:sectPr>'
           f'</w:body></w:document>')
    ct = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
          '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
          '<Default Extension="xml" ContentType="application/xml"/>'
          '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-'
          'officedocument.wordprocessingml.document.main+xml"/></Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
            'relationships/officeDocument" Target="word/document.xml"/></Relationships>')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", ct)
        z.writestr("_rels/.rels", rels)
        z.writestr("word/document.xml", doc)


def main():
    full = {o[0]: load_rav(o[3], o[5]) for o in OBS if o[0] in COLS + ("rbdA",)}
    win = {k: v[:, -N_LAST:] for k, v in full.items()}
    up = win["rbdA"].mean(1) > RBD_UP_CUTOFF

    rows = []
    for k, i in enumerate(SPIKES):
        cells = []
        for c in COLS:
            m, sd, f0 = win[c][k].mean(), win[c][k].std(ddof=1), full[c][k, 0]
            if c == "ntd":
                cells += [f"{f0:.0f}", f"{m:.0f} ± {sd:.0f}"]
            else:
                cells += [f"{f0:.1f}", f"{m:.1f} ± {sd:.1f}"]
        rows.append([str(i)] + cells + ["up" if up[k] else "down"])

    # ensemble row: mean ± SD of the 29 per-spike means
    ens = []
    for c in COLS:
        m, f0 = win[c].mean(1), full[c][:, 0]
        if c == "ntd":
            ens += [f"{f0.mean():.0f} ± {f0.std(ddof=1):.0f}",
                    f"{m.mean():.0f} ± {m.std(ddof=1):.0f}"]
        else:
            ens += [f"{f0.mean():.1f} ± {f0.std(ddof=1):.1f}",
                    f"{m.mean():.1f} ± {m.std(ddof=1):.1f}"]
    ens_row = ["All 29"] + ens + [f"{int(up.sum())} up / {int((~up).sum())} down"]

    header = ["Spike"] + list(HEAD)

    with open(os.path.join(HERE, "SI_table_per_spike.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
        w.writerow(ens_row)

    with open(os.path.join(HERE, "SI_table_per_spike.md"), "w") as f:
        f.write("| " + " | ".join(header) + " |\n")
        f.write("|" + "|".join(["---"] * len(header)) + "|\n")
        for r in rows:
            f.write("| " + " | ".join(r) + " |\n")
        f.write("| **" + "** | **".join(ens_row) + "** |\n")

    with open(os.path.join(HERE, "SI_table_per_spike.tex"), "w") as f:
        f.write("\\begin{tabular}{l" + "r" * (len(header) - 1) + "}\n\\toprule\n")
        f.write(" & ".join(h.replace("°", "$^\\circ$").replace("Å²", "\\AA$^2$")
                           for h in header) + " \\\\\n\\midrule\n")
        for r in rows:
            f.write(" & ".join(c.replace("±", "$\\pm$") for c in r) + " \\\\\n")
        f.write("\\midrule\n")
        f.write(" & ".join(c.replace("±", "$\\pm$") for c in ens_row) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")

    with open(os.path.join(HERE, "SI_table_per_spike.html"), "w") as f:
        f.write("<html><head><meta charset='utf-8'><style>"
                "table{border-collapse:collapse;font-family:Calibri,sans-serif;font-size:10pt}"
                "th,td{border:0;padding:3px 10px;text-align:right}"
                "th:first-child,td:first-child{text-align:left}"
                "thead tr{border-top:1px solid #000;border-bottom:1px solid #000}"
                "tbody tr:last-child{border-top:1px solid #000;font-weight:bold}"
                "table{border-bottom:1px solid #000}"
                "</style></head><body><table><thead><tr>")
        f.write("".join(f"<th>{h}</th>" for h in header) + "</tr></thead><tbody>")
        for r in rows:
            f.write("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>")
        f.write("<tr>" + "".join(f"<td>{c}</td>" for c in ens_row) + "</tr>")
        f.write("</tbody></table></body></html>\n")

    write_docx(os.path.join(HERE, "SI_table_per_spike.docx"), header, rows, ens_row)

    print(f"wrote SI_table_per_spike.{{csv,md,tex,html,docx}} ({len(rows)} spikes)")


if __name__ == "__main__":
    main()
