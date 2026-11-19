"""Build the editable article, supporting information and typeset reading copies."""
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
MARKDOWN = r"markdown+tex_math_single_backslash"
IMAGE_AND_CAPTION = re.compile(r"^!\[([^\]]*)\]\((figures/[^)]+)\)\n\n(\*\*(?:Supplementary )?Figure [^\n]*)", re.M)


def run(*args, **kwargs):
    return subprocess.run(args, cwd=HERE, check=True, **kwargs)


def convert_caption(text):
    return run("pandoc", "--from", MARKDOWN, "--to", "latex", input=text,
               text=True, capture_output=True).stdout.strip()


def pdf_source(text, wide):
    def replace(match):
        path = Path(match[2]).with_suffix(".pdf").as_posix()
        env = "figure*" if wide else "figure"
        return (
            "\n```{=latex}\n"
            rf"\begin{{{env}}}[!htbp]" + "\n\\centering\n"
            rf"\includegraphics[width=\textwidth]{{{path}}}" + "\n"
            r"\par\smallskip\begin{minipage}{\textwidth}\small" + "\n"
            + convert_caption(match[3]) + "\n"
            r"\end{minipage}" + "\n" + rf"\end{{{env}}}" + "\n```\n"
        )
    return IMAGE_AND_CAPTION.sub(replace, text)


def normalize_typesetting(source):
    # Use typographic super/subscripts with the chosen text font and DOCX.
    for characters, replacements, marker in [
        ("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789−+", "^"),
        ("₀₁₂₃₄₅₆₇₈₉ₚₛ", "0123456789ps", "~"),
    ]:
        table = str.maketrans(characters, replacements)
        source = re.sub(f"[{characters}]+",
                        lambda match: marker + match[0].translate(table) + marker,
                        source)
    return re.sub(r"\{\\rm\s+([^{}]+)\}", r"\\mathrm{\1}", source).replace("∈", r" \(\in\) ")


def main():
    run(sys.executable, str(PAPER / "check_callout_order.py"))
    draft = normalize_typesetting((PAPER / "MANUSCRIPT_WORKING_DRAFT_EN.md").read_text())
    article, rest = draft.split("## Supplementary Results", 1)
    supplemental, refs = rest.split("## References", 1)
    title, remainder = article.split("\n", 1)
    abstract, text = remainder.split("## Introduction", 1)
    abstract = abstract.split("## Abstract", 1)[1].strip()
    significance = (HERE / "SIGNIFICANCE.md").read_text().strip()
    assert len(abstract.split()) <= 300 and len(significance.split()) <= 120

    byline = (HERE / "AUTHOR_FRONT_MATTER.md").read_text().strip()
    end = """## Funding

Renato Socodato holds an employment contract funded by the Fundação para a Ciência e a Tecnologia (FCT), Portugal (2023.06380.CEECIND).

## AI disclosure

AI-assisted tools were used for code review and manuscript proofreading.

## Data and code availability

The article reports analytical results and simulations. Previously published measurements used for comparison are cited where they are discussed. Code and figure reproduction materials are maintained in a private repository pending release.

"""
    article_body = "## Introduction" + text + end + "## References" + refs
    supplementary_refs = re.sub(r"(?m)^26\. CD11b-enriched mouse brain proteome[^\n]*\n?", "", refs)
    supplement = (
        "# Supporting information\n\n" + title.removeprefix("# ") + "\n\n"
        "## Supplementary Results" + supplemental + "\n## References used in the article and supplement" + supplementary_refs
    )
    article_body = article_body.replace("submission_biophysj/figures/", "figures/")
    supplement = supplement.replace("submission_biophysj/figures/", "figures/")
    article = (
        title + "\n\n" + byline + "\n\n## Abstract\n\n" + abstract
        + "\n\n## Significance\n\n" + significance + "\n\n" + article_body
    )
    (HERE / "manuscript.md").write_text(article)
    (HERE / "supporting_information.md").write_text(supplement)
    if (HERE / "COVER_LETTER.md").exists():
        run("pandoc", "COVER_LETTER.md", "--from", MARKDOWN, "-o", "cover_letter.docx")
    for stem, content, wide in [
        ("manuscript", article, True), ("supporting_information", supplement, False)
    ]:
        run("pandoc", f"{stem}.md", "--from", MARKDOWN, "--resource-path=.", "-o", f"{stem}.docx")
        if wide:
            heading = convert_caption(title.removeprefix("# "))
            abs_tex = convert_caption(abstract)
            sig_tex = convert_caption(significance)
            front = (
                "```{=latex}\n\\twocolumn[{\\begin{minipage}{\\textwidth}\n"
                "{\\sffamily\\fontsize{9}{11}\\selectfont Regular Article}\\par\\vspace{2mm}\n"
                "{\\raggedright\\sffamily\\bfseries\\fontsize{20}{22}\\selectfont " + heading + "\\par}\\vspace{2mm}\n"
                "{\\sffamily\\fontsize{10}{12}\\selectfont " + convert_caption(byline) + "}\\par\\vspace{3mm}\n"
                "{\\sffamily\\bfseries\\fontsize{9}{11}\\selectfont ABSTRACT}\\quad " + abs_tex + "\\par\\medskip\n"
                "{\\sffamily\\bfseries\\fontsize{9}{11}\\selectfont SIGNIFICANCE}\\quad " + sig_tex
                + "\\par\\vspace{4mm}\\end{minipage}}]\n```\n\n"
            )
            latex = front + pdf_source(article_body, wide=True)
            latex = latex.replace('## References', '```{=latex}\n\\end{multicols}\n\\clearpage\n```\n\n' + '## References\n\n' + '```{=latex}\n\\begin{multicols}{2}\\small\n```\n\n', 1)
            latex += '\n\n' + '```{=latex}\n\\end{multicols}\n```\n\n'
            latex = latex.replace("## Discussion", '```{=latex}\n\\onecolumn\n\\begin{multicols}{2}\n```\n\n' + "## Discussion", 1)
        else:
            latex = pdf_source(supplement, wide=False)
            latex = latex.replace("## References used", "```{=latex}\n\\small\n```\n\n## References used", 1)
        render_source = HERE / f"{stem}_typeset.md"
        render_source.write_text(latex)
        run("pandoc", render_source.name, "--from", MARKDOWN, "--template", "bj_reading.latex",
            "--pdf-engine=xelatex", *(["-V", "onecolumn=true"] if not wide else []), "-o", f"{stem}.pdf")


if __name__ == "__main__":
    main()
