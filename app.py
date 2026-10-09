import ctypes
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from htmltools import HTMLDependency
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shiny.render.renderer import Renderer
from shiny.types import SafeException
from scipy.stats import linregress
import numpy as np

genes = pd.read_csv("DETECTED_GENES.csv", header=None, index_col=False)

ages = ["Young", "Adult", "Old"]
sexes = ["Female", "Male"]
lines = ["Microglia", "Neurons", "Astrocytes"]
comp_order = ["modCG", "mCG", "hmCG"]

def sort_key(g):
    starts_with_digit = g[0].isdigit()
    ends_with_rik = g.lower().endswith('rik')
    return (ends_with_rik, starts_with_digit, g.lower())

sorted_genes = sorted(genes[0], key=sort_key)
gene_set = set(sorted_genes)

color_map = {
    'Young': "#74a5ce",
    'Adult': "#afc75b",
    'Old': "#b14646",
    'Female': "#b33ba9",
    'Male': "#4c3ec7",
    'Microglia': "#bb58da",
    'Neurons': "#53410F",
    'Astrocytes': "#4d26bb",
    'Young Female': "#f088e7",
    'Young Male': "#4ecbeb",
    'Adult Female': "#c65cdb",
    'Adult Male': "#6d93e6",
    'Old Female': "#eb25eb",
    'Old Male': "#4224ee",
    'Young Microglia': "#d894f3",
    'Young Neurons': "#ecd279",
    'Young Astrocytes': "#9694f3",
    'Adult Microglia': "#bc58d4",
    'Adult Neurons': "#ebb238",
    'Adult Astrocytes': "#6b58d4",
    'Old Microglia': "#cd0ff3",
    'Old Neurons': "#b37708",
    'Old Astrocytes': "#130ff3",
    'Female Microglia': "#dc24ec",
    'Male Microglia': "#a351f0",
    'Female Neurons': "#e9a617",
    'Male Neurons': "#e45311",
    'Female Astrocytes': "#7708f5",
    'Male Astrocytes': "#0e49eb",
    'Young Female Microglia': "#e5baee",
    'Young Male Microglia': "#89d8a7",
    'Adult Female Microglia': "#a550c7",
    'Adult Male Microglia': "#57A063",
    'Old Female Microglia': "#4c0d5f",
    'Old Male Microglia': "#043608",
    'Young Female Neurons': "#f0d68e",
    'Young Male Neurons': "#8ad6c6",
    'Adult Female Neurons': "#926C3B",
    'Adult Male Neurons': "#319988",
    'Old Male Neurons': "#094250",
    'Old Female Neurons': "#5f3404",
    'Young Female Astrocytes': "#eca57c",
    'Young Male Astrocytes': "#97b8e2",
    'Adult Female Astrocytes': "#eb8d40",
    'Adult Male Astrocytes': "#49b3e4",
    'Old Female Astrocytes': "#e44040",
    'Old Male Astrocytes': "#4839c9"
}

# "Grouping" choice -> factors that are filtered on and combined into the group label
groupings = {
    "1": ["AGE"],
    "2": ["SEX"],
    "3": ["LINE"],
    "4": ["AGE", "SEX"],
    "5": ["AGE", "LINE"],
    "6": ["SEX", "LINE"],
    "7": ["AGE", "SEX", "LINE"],
}

def load(path, id_col):
    data = pq.read_table(path).to_pandas()
    data[id_col] = data[id_col].astype(str)
    data["AGE"] = pd.Categorical(data["AGE"], categories=ages, ordered=True)
    data["SEX"] = pd.Categorical(data["SEX"], categories=sexes, ordered=True)
    data["LINE"] = pd.Categorical(data["LINE"], categories=lines, ordered=True)
    sort_cols = ['LINE', 'AGE', 'SEX']
    if "COMP" in data.columns:
        data["COMP"] = pd.Categorical(data["COMP"], categories=comp_order, ordered=True)
        sort_cols.append('COMP')
    return data.sort_values(sort_cols)

# Load all data once at startup. Reading parquet per request is slow and blocks every
# session, since Shiny runs this code on a single event loop.
expr_data = load('ALL_TPM_LABELED_FILTERED.parquet', "geneid")
body_data = load('ALL_GENE_BODY_PER_SAMPLE_T_v2.parquet', "gene")
tss_data = load('ALL_TSS_PER_SAMPLE_T_v2.parquet', "gene")

# Parsing the very wide parquet footers leaves ~350 MB of freed heap that glibc keeps
# unless asked to give it back to the OS (Linux only; harmless to skip elsewhere).
try:
    ctypes.CDLL("libc.so.6").malloc_trim(0)
except (OSError, AttributeError):
    pass

def merge_corr(mod_data, gene):
    # Pair each sample's expression (<gene>_x) with its modification levels (<gene>_y)
    data = pd.merge(
        expr_data[["geneid", "AGE", "SEX", "LINE", gene]],
        mod_data[["gene", "AGE", "SEX", "LINE", "COMP", gene]],
        left_on = ['geneid', 'AGE', 'SEX', 'LINE'], right_on = ['gene', 'AGE', 'SEX', 'LINE'])
    return data.sort_values(['LINE', 'AGE', 'SEX', 'COMP'])

def apply_grouping(data, factors, selected, value_cols):
    """Filter rows to the selected levels of each factor and add a GROUP label.

    `selected` maps each factor to a function returning its selected levels, so only the
    inputs for the active grouping become reactive dependencies.
    """
    keep = [f for f in ("AGE", "SEX", "LINE") if f in factors or f == "LINE"]
    # Select columns before filtering rows, so the wide tables aren't copied
    outputData = data[keep + value_cols]
    for factor in factors:
        outputData = outputData.loc[outputData[factor].isin(selected[factor]())]
    outputData = outputData.copy()
    group = outputData[factors[0]].astype(str)
    for factor in factors[1:]:
        group = group + " " + outputData[factor].astype(str)
    outputData["GROUP"] = group
    return outputData

def theme(dark):
    """Return (plot background, font color, paper background)."""
    if dark:
        return "#B9B9B9", "#ffffff", "#252525"
    return "#e4e4e4", "#000000", "#ffffff"

# Figures are built as plain dicts and drawn by plotly.js in the browser. Building
# plotly.graph_objects figures and sending them as widgets took ~0.6 s per plot, all of
# it blocking every session; plain dicts take a few milliseconds.
def base_layout(fig):
    """Default layout of `fig` (template and subplot positions included) as a dict.

    Built once with plotly so the plots look the same as plotly-built figures. The
    template margins match what shinywidgets used.
    """
    fig.layout.template.layout.margin = dict(l=16, t=32, r=16, b=16)
    return fig.to_plotly_json()["layout"]

single_layout = base_layout(go.Figure())
three_panel_layout = base_layout(make_subplots(rows=1, cols=3))

def build_layout(base, **updates):
    """Copy `base`, merging dict values one level deep (e.g. adding an axis title)."""
    layout = dict(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(layout.get(key), dict):
            layout[key] = {**layout[key], **value}
        else:
            layout[key] = value
    return layout

def panel_axes(position):
    """Axis references for a trace in panel 1-3 of a three-panel figure."""
    suffix = "" if position == 1 else str(position)
    return {"xaxis": "x" + suffix, "yaxis": "y" + suffix}

def json_values(values):
    """Floats as a JSON-safe list; NaN becomes None, which plotly.js treats as missing."""
    return [None if v != v else v for v in np.asarray(values, dtype=float).tolist()]

def title_font(color):
    return dict(size = 24, textcase = "upper", weight = "bold", color = color)

def figure(data, layout):
    # Unless the height is fixed, fill the output container and follow its size
    # (as shinywidgets did)
    config = {} if "height" in layout else {"responsive": True}
    return {"data": data, "layout": layout, "config": config}

plotly_dependencies = [
    HTMLDependency(
        "plotly", plotly.offline.get_plotlyjs_version(),
        source={"subdir": str(Path(plotly.__file__).parent / "package_data")},
        script={"src": "plotly.min.js"}, all_files=False),
    HTMLDependency(
        "plotly-output", "1.0",
        source={"subdir": str(Path(__file__).parent / "www")},
        script={"src": "plotly_output.js"}, stylesheet={"href": "plotly_output.css"}, all_files=False),
]

def output_plotly(id):
    """Placeholder for a @render_plotly figure; fills its container like output_widget()."""
    return ui.fill.as_fillable_container(ui.fill.as_fill_item(
        ui.div(*plotly_dependencies, id=id, class_="plotly-output")))

class render_plotly(Renderer[dict]):
    """Send a figure dict ({data, layout, config}) to the browser (see www/plotly_output.js)."""
    async def transform(self, value: dict) -> dict:
        return value

app_ui = ui.page_sidebar(
    ui.sidebar(
        # Only the default gene is written into the page; the full list (~25,000 genes)
        # is searched on the server as the visitor types (see update_selectize in
        # server()). Writing every gene into the page froze the browser on load.
        ui.input_selectize(
            "gene",
            "Gene",
            ["Cx3cr1"],
            selected="Cx3cr1"
        ),
        ui.input_select(
            "filter",
            "Grouping",
            {
                1: "Age",
                2: "Sex",
                3: "Cell Type",
                4: "Age and Sex",
                5: "Age and Cell Type",
                6: "Sex and Cell Type",
                7: "Age, Sex, and Cell Type"
            },
            selected=7
        ),
        ui.panel_conditional(
            "input.filter == 1 || input.filter == 4 || input.filter == 5 || input.filter == 7",
            ui.input_checkbox_group(
                "age",
                "Filter by Age",
                {
                    "Young": "Young (3 mo)",
                    "Adult": "Adult (12 mo)",
                    "Old": "Old (24 mo)"
                },
                selected=ages
            ),
        ),
        ui.panel_conditional(
            "input.filter == 2 || input.filter == 4 || input.filter == 6 || input.filter == 7",
            ui.input_checkbox_group(
                "sex",
                "Filter by Sex",
                {
                    "Female": "Female",
                    "Male": "Male"
                },
                selected=sexes
            ),
        ),
        ui.panel_conditional(
            "input.filter == 3 || input.filter == 5 || input.filter == 6 || input.filter == 7",
            ui.input_checkbox_group(
                "line",
                "Filter by Cell Type",
                ["Astrocytes", "Neurons", "Microglia"],
                selected=lines
            ),
        ),
        ui.download_button("download_expr", "Download Expression Data"),
        ui.download_button("download_gene", "Download Gene Body Modification Data"),
        ui.download_button("download_tss", "Download Promoter Modification Data"),
        ui.input_action_button("toggle_dark", "Toggle Dark Mode")
    ),
    ui.page_navbar(
        ui.nav_panel(
            "Gene Expression",
            ui.layout_columns(
                output_plotly("expression_plot")
            )
        ),
        ui.nav_panel(
            "Gene DNA Modifications",
            ui.layout_columns(
                output_plotly("tss_plot")
            ),
            ui.layout_columns(
                output_plotly("gene_body_plot")
            )
        ),
        ui.nav_panel(
            "Correlation Plots",
            ui.layout_columns(
                output_plotly("tss_corr_plot")
            ),
            ui.layout_columns(
                output_plotly("gene_corr_plot")
            )

        )
    ),
    ui.hr(),
    ui.tags.footer("The promoter region is defined as 500 bases upstream to 100 bases downstream of the TSS."),
    ui.tags.footer("The gene body is defined as the TSS to the TES."),
    ui.tags.footer("The study associated with this data can be found here: .")
)

def server(input: Inputs, output: Outputs, session: Session):
    # Per-session, so one visitor's toggle doesn't affect anyone else
    dark = reactive.value(False)

    ui.update_selectize("gene", choices=sorted_genes, selected="Cx3cr1", server=True)

    selected = {"AGE": input.age, "SEX": input.sex, "LINE": input.line}

    @render.download_button(filename="gene_expression_data.csv")
    def download_expr():
        outData = filtered_expr()
        yield outData.to_csv(index=False)

    @render.download_button(filename="gene_body_methylation_data.csv")
    def download_gene():
        outData = filtered_body()
        yield outData.to_csv(index=False)

    @render.download_button(filename="promoter_methylation_data.csv")
    def download_tss():
        outData = filtered_tss()
        yield outData.to_csv(index=False)

    @reactive.effect
    @reactive.event(input.toggle_dark)
    def _():
        dark.set(not dark())
        ui.update_dark_mode("dark" if dark() else "light")

    # Browsers can send any value for an input, not just the offered choices,
    # so validate before using them to index the data.
    @reactive.calc
    def gene() -> str:
        g = input.gene()
        req(isinstance(g, str) and g in gene_set)
        return g

    @reactive.calc
    def factors() -> list[str]:
        f = input.filter()
        req(isinstance(f, str) and f in groupings)
        return groupings[f]

    def require_gene_in(data, description):
        if gene() not in data.columns:
            raise SafeException(f"No {description} data is available for {gene()}.")

    @reactive.calc
    def filtered_expr() -> pd.DataFrame:
        require_gene_in(expr_data, "expression")
        return apply_grouping(expr_data, factors(), selected, [gene()])

    @render_plotly
    def expression_plot():
        data = filtered_expr()
        bgcolor, fontcolor, papercolor = theme(dark())
        traces = []
        for group in data['GROUP'].unique():
            traces.append(dict(type = "box", y = json_values(data[data['GROUP'] == group][gene()]),
                boxpoints = 'all', jitter = 0.5, marker = dict(line = dict(width=1), color = color_map[group]),
                line = dict(width=2), pointpos = 0, name = group))
        layout = build_layout(single_layout,
            title = dict(text = gene(), font = title_font(fontcolor)),
            font = dict(color = fontcolor),
            yaxis = dict(title = dict(text = "TPM")),
            xaxis = dict(title = dict(text = "")),
            showlegend = False,
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin = dict(t=100)
        )
        return figure(traces, layout)

    @reactive.calc
    def filtered_body() -> pd.DataFrame:
        require_gene_in(body_data, "gene body modification")
        return apply_grouping(body_data, factors(), selected, ["COMP", gene()])

    @reactive.calc
    def filtered_tss() -> pd.DataFrame:
        require_gene_in(tss_data, "promoter modification")
        return apply_grouping(tss_data, factors(), selected, ["COMP", gene()])

    def modification_plot(mod_data, title, clip_lower):
        position=1
        bgcolor, fontcolor, papercolor = theme(dark())
        traces = []
        for comp in mod_data['COMP'].unique():
            comp_data = mod_data[mod_data['COMP'] == comp]
            for group in comp_data['GROUP'].unique():
                traces.append(dict(type = "box", y = json_values(comp_data[comp_data['GROUP'] == group][gene()].clip(lower=clip_lower)),
                    boxpoints = 'all', jitter = 0.5, marker = dict(line = dict(width=1), color = color_map[group]),
                    line = dict(width=1), pointpos = 0, name = group, **panel_axes(position)))
            position += 1
        layout = build_layout(three_panel_layout,
            title = dict(text = title + gene(), font = title_font(fontcolor)),
            showlegend = False,
            font = dict(color = fontcolor),
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin = dict(t=100),
            yaxis = dict(title = dict(text = "5modCG (%)")),
            yaxis2 = dict(title = dict(text = "5mCG (%)")),
            yaxis3 = dict(title = dict(text = "5hmCG (%)"))
        )
        return figure(traces, layout)

    @render_plotly
    def gene_body_plot():
        return modification_plot(filtered_body(), "Gene Body Methylation: ", clip_lower=0)

    @render_plotly
    def tss_plot():
        return modification_plot(filtered_tss(), "Promoter Methylation: ", clip_lower=10)

    @reactive.calc
    def filtered_gene_corr() -> pd.DataFrame:
        require_gene_in(expr_data, "expression")
        require_gene_in(body_data, "gene body modification")
        data = merge_corr(body_data, gene())
        return apply_grouping(data, factors(), selected, [f"{gene()}_x", f"{gene()}_y", "COMP"])

    @reactive.calc
    def filtered_tss_corr() -> pd.DataFrame:
        require_gene_in(expr_data, "expression")
        require_gene_in(tss_data, "promoter modification")
        data = merge_corr(tss_data, gene())
        return apply_grouping(data, factors(), selected, [f"{gene()}_x", f"{gene()}_y", "COMP"])

    def correlation_plot(corr_data, title):
        position = 1
        bgcolor, fontcolor, papercolor = theme(dark())
        traces = []

        legendKey = [False, False, True]

        for comp in corr_data['COMP'].unique():
            comp_data = corr_data[corr_data['COMP'] == comp]

            # ---- Scatter points for each group ----
            for group in comp_data['GROUP'].unique():
                subset = comp_data[comp_data['GROUP'] == group]

                x = subset[f"{gene()}_x"]
                y = subset[f"{gene()}_y"]

                # Remove NaNs and ensure x>0 for log scale
                mask = x.notna() & y.notna() & (x > 0)
                x = x[mask]
                y = y[mask]

                # Floor y values at 0
                y = np.maximum(y, 0)

                traces.append(dict(
                    type="scatter",
                    x=json_values(x),
                    y=json_values(y),
                    name=group,
                    marker=dict(color=color_map[group]),
                    mode="markers",
                    legendgroup="test",
                    showlegend=legendKey[position-1],
                    **panel_axes(position)
                ))

            # ---- One regression line across all groups ----
            x_all = comp_data[f"{gene()}_x"]
            y_all = comp_data[f"{gene()}_y"]

            mask_all = x_all.notna() & y_all.notna() & (x_all > 0)
            x_all = x_all[mask_all]
            y_all = np.maximum(y_all[mask_all], 0)

            if len(x_all) > 1:
                x_log = np.log10(x_all)
                result = linregress(x_log, y_all)

                x_sorted = np.sort(x_all)
                y_fit = result.slope * np.log10(x_sorted) + result.intercept
                y_fit = np.maximum(y_fit, 0)  # keep regression line >=0

                traces.append(dict(
                    type="scatter",
                    x=json_values(x_sorted),
                    y=json_values(y_fit),
                    mode="lines",
                    line=dict(color="black", width=2),
                    name="Trendline",
                    showlegend=False,
                    **panel_axes(position)
                ))

            position += 1

        log_tpm_axis = dict(type="log", title=dict(text="TPM"))
        layout = build_layout(three_panel_layout,
            title=dict(text=title + gene(), font=title_font(fontcolor)),
            showlegend=True,
            font=dict(color=fontcolor),
            paper_bgcolor=papercolor,
            plot_bgcolor=bgcolor,
            margin=dict(t=100),
            height=600,
            yaxis=dict(title=dict(text="5modCG (%)")),
            yaxis2=dict(title=dict(text="5mCG (%)")),
            yaxis3=dict(title=dict(text="5hmCG (%)")),
            xaxis=log_tpm_axis, xaxis2=log_tpm_axis, xaxis3=log_tpm_axis
        )

        return figure(traces, layout)

    @render_plotly
    def gene_corr_plot():
        return correlation_plot(filtered_gene_corr(), "Gene Body Correlation: ")

    @render_plotly
    def tss_corr_plot():
        return correlation_plot(filtered_tss_corr(), "Promoter Correlation: ")

app = App(app_ui, server)
# Show a generic message instead of raw exception text in the browser
# (SafeException messages, like "no data for this gene", are still shown).
app.sanitize_errors = True
