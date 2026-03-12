import pyarrow.parquet as pq
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shinywidgets import output_widget, render_widget
from scipy.stats import linregress
import numpy as np

genes = pd.read_csv("DETECTED_GENES.csv", header=None, index_col=False)

mode="light"

ages = ["Young", "Adult", "Old"]
sexes = ["Female", "Male"]
lines = ["Microglia", "Neurons", "Astrocytes"]
comp_order = ["modCG", "mCG", "hmCG"]

def sort_key(g):
    starts_with_digit = g[0].isdigit()
    ends_with_rik = g.lower().endswith('rik')
    return (ends_with_rik, starts_with_digit, g.lower())

sorted_genes = sorted(genes[0], key=sort_key)

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

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_selectize(
            "gene", 
            "Gene", 
            sorted_genes,
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
        ui.download_button("download_gene", "Download Gene Body Modificaiton Data"),
        ui.download_button("download_tss", "Download Promoter Modificaiton Data"),
        ui.input_action_button("toggle_dark", "Toggle Dark Mode")
    ),
    ui.page_navbar(
        ui.nav_panel(
            "Gene Expression",
            ui.layout_columns(
                output_widget("expression_plot")
            )
        ),
        ui.nav_panel(
            "Gene DNA Modifications",
            ui.layout_columns(
                output_widget("gene_body_plot")
            ),
            ui.layout_columns(
                output_widget("tss_plot")
            )
        ),
        ui.nav_panel(
            "Correlation Plots",
            ui.layout_columns(
                output_widget("gene_corr_plot")
            ),
            ui.layout_columns(
                output_widget("tss_corr_plot")
            )
        )
    )
)

def server(input: Inputs, output: Outputs, session: Session):
    @render.download(filename="gene_expression_data.csv")
    def download_expr():
        outData = filtered_expr()
        yield outData.to_csv(index=False)
    
    @render.download(filename="gene_body_methylation_data.csv")
    def download_gene():
        outData = filtered_body()
        yield outData.to_csv(index=False)

    @render.download(filename="promoter_methylation_data.csv")
    def download_tss():
        outData = filtered_tss()
        yield outData.to_csv(index=False)
    
    @reactive.effect
    @reactive.event(input.toggle_dark)
    def _():
        global mode
        if mode == "light":
            ui.update_dark_mode("dark")
            mode = "dark"
        else:
            ui.update_dark_mode("light")
            mode = "light"
    
    @reactive.Calc
    def filtered_expr() -> pd.DataFrame:
        data = pq.read_table('ALL_RPKM_LABELED_FILTERED.parquet', columns=["AGE", "SEX", "LINE", input.gene()]).to_pandas()
        data["AGE"] = pd.Categorical(data["AGE"], categories=ages, ordered=True)
        data["SEX"] = pd.Categorical(data["SEX"], categories=sexes, ordered=True)
        data["LINE"] = pd.Categorical(data["LINE"], categories=lines, ordered=True)
        sorted_data = data.sort_values(['LINE', 'AGE', 'SEX'])
        match input.filter():
            case "1":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "LINE", input.gene())]
                outputData["GROUP"] = outputData["AGE"]
                return outputData
            case "2":
                outputData = sorted_data.loc[sorted_data["SEX"].isin(input.sex()), ("SEX", "LINE", input.gene())]
                outputData["GROUP"] = outputData["SEX"]
                return outputData
            case "3":
                outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("LINE", input.gene())]
                outputData["GROUP"] = outputData["LINE"]
                return outputData
            case "4":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
                return outputData
            case "5":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "LINE", input.gene())]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "6":
                outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("SEX", "LINE", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "7":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
                return outputData

    @render_widget
    def expression_plot():
        if input.gene() not in sorted_genes:
            return
        fig = go.Figure()
        data = filtered_expr()
        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"
        for group in data['GROUP'].unique():
            fig.add_trace(go.Box(y = data[data['GROUP'] == group][input.gene()],
                boxpoints = 'all', jitter = 0.5, marker_line_width=1, line = dict(width=2),
                pointpos = 0, name = group, marker_color=color_map[group]))
        fig.update_layout(
            title=input.gene(),
            title_font = dict(
                size = 24,
                textcase = "upper",
                weight = "bold",
                color = fontcolor
            ),
            font = dict(
                color = fontcolor
            ),
            yaxis_title = "RPKM",
            xaxis_title = "",
            showlegend = False,
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin=dict(t=100)
        )
        return fig
    
    @reactive.Calc
    def filtered_body() -> pd.DataFrame:
        gene_corr = pq.read_table('ALL_GENE_BODY_PER_SAMPLE_TRANSPOSED_FINAL.parquet', columns=["AGE", "SEX", "LINE", "COMP", input.gene()]).to_pandas()
        gene_corr["AGE"] = pd.Categorical(gene_corr["AGE"], categories=ages, ordered=True)
        gene_corr["SEX"] = pd.Categorical(gene_corr["SEX"], categories=sexes, ordered=True)
        gene_corr["LINE"] = pd.Categorical(gene_corr["LINE"], categories=lines, ordered=True)
        gene_corr["COMP"] = pd.Categorical(gene_corr["COMP"], categories=comp_order, ordered=True)
        sorted_body = gene_corr.sort_values(['LINE', 'AGE', 'SEX', 'COMP'])
        match input.filter():
            case "1":
                outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
                outputData["GROUP"] = outputData["AGE"]
                return outputData
            case "2":
                outputData = sorted_body.loc[sorted_body["SEX"].isin(input.sex()), ("SEX", "LINE", "COMP", input.gene())]
                outputData["GROUP"] = outputData["SEX"]
                return outputData
            case "3":
                outputData = sorted_body.loc[sorted_body["LINE"].isin(input.line()), ("LINE", "COMP", input.gene())]
                outputData["GROUP"] = outputData["LINE"]
                return outputData
            case "4":
                outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
                return outputData
            case "5":
                outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "6":
                outputData = sorted_body.loc[sorted_body["LINE"].isin(input.line()), ("SEX", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "7":
                outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
                return outputData

    @render_widget
    def gene_body_plot():
        if input.gene() not in sorted_genes:
            return
        fig = make_subplots(rows=1, cols=3)
        body_data = filtered_body()
        position=1
        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"
        for comp in body_data['COMP'].unique():
            comp_data = body_data[body_data['COMP'] == comp]
            for group in comp_data['GROUP'].unique():
                fig.add_trace(go.Box(y = comp_data[comp_data['GROUP'] == group][input.gene()],
                    boxpoints = 'all', jitter = 0.5, marker_line_width=1, line = dict(width=1),
                    pointpos = 0, name = group, marker_color=color_map[group]),
                    row=1, col=position)
            position += 1
        fig.update_layout(
            title = "Gene Body Methylation: " + input.gene(),
            showlegend = False,
            title_font = dict(
                size = 24,
                textcase = "upper",
                weight = "bold",
                color = fontcolor
            ),
            font = dict(
                color = fontcolor
            ),
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin=dict(t=100)
        )
        # range=[0, 100]
        fig.update_yaxes(title_text="5modCG (%)", row=1, col=1)
        fig.update_yaxes(title_text="5mCG (%)", row=1, col=2)
        fig.update_yaxes(title_text="5hmCG (%)", row=1, col=3)
        return fig

    @reactive.Calc
    def filtered_tss() -> pd.DataFrame:
        tss_corr = pq.read_table('ALL_TSS_PER_SAMPLE_TRANSPOSED.parquet', columns=["AGE", "SEX", "LINE", "COMP", input.gene()]).to_pandas()
        tss_corr["AGE"] = pd.Categorical(tss_corr["AGE"], categories=ages, ordered=True)
        tss_corr["SEX"] = pd.Categorical(tss_corr["SEX"], categories=sexes, ordered=True)
        tss_corr["LINE"] = pd.Categorical(tss_corr["LINE"], categories=lines, ordered=True)
        tss_corr["COMP"] = pd.Categorical(tss_corr["COMP"], categories=comp_order, ordered=True)
        sorted_tss = tss_corr.sort_values(['LINE', 'AGE', 'SEX', 'COMP'])
        match input.filter():
            case "1":
                outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
                outputData["GROUP"] = outputData["AGE"]
                return outputData
            case "2":
                outputData = sorted_tss.loc[sorted_tss["SEX"].isin(input.sex()), ("SEX", "LINE", "COMP", input.gene())]
                outputData["GROUP"] = outputData["SEX"]
                return outputData
            case "3":
                outputData = sorted_tss.loc[sorted_tss["LINE"].isin(input.line()), ("LINE", "COMP", input.gene())]
                outputData["GROUP"] = outputData["LINE"]
                return outputData
            case "4":
                outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
                return outputData
            case "5":
                outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "6":
                outputData = sorted_tss.loc[sorted_tss["LINE"].isin(input.line()), ("SEX", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "7":
                outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
                return outputData

    @render_widget
    def tss_plot():
        if input.gene() not in sorted_genes:
            return
        fig = make_subplots(rows=1, cols=3)
        tss_data = filtered_tss()
        position=1
        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"
        for comp in tss_data['COMP'].unique():
            comp_data = tss_data[tss_data['COMP'] == comp]
            for group in comp_data['GROUP'].unique():
                fig.add_trace(go.Box(y = comp_data[comp_data['GROUP'] == group][input.gene()],
                    boxpoints = 'all', jitter = 0.5, marker_line_width=1, line = dict(width=1),
                    pointpos = 0, name = group, marker_color=color_map[group]),
                    row=1, col=position)
            position += 1
        fig.update_layout(
            title = "Promoter Methylation: " + input.gene(),
            showlegend = False,
            title_font = dict(
                size = 24,
                textcase = "upper",
                weight = "bold",
                color = fontcolor
            ),
            font = dict(
                color = fontcolor
            ),
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin=dict(t=100)
        )
        # range=[0, 100]
        fig.update_yaxes(title_text="5modCG (%)", row=1, col=1)
        fig.update_yaxes(title_text="5mCG (%)", row=1, col=2)
        fig.update_yaxes(title_text="5hmCG (%)", row=1, col=3)
        return fig
    
    @reactive.Calc
    def filtered_gene_corr() -> pd.DataFrame:
        rpkm_corr = pq.read_table('ALL_RPKM_DATA_FILTERED_T_v2.parquet', columns=["gene", "AGE", "SEX", "LINE", input.gene()]).to_pandas()
        gene_corr = pq.read_table('ALL_GENE_BODY_PER_SAMPLE_T_v2.parquet', columns=["gene", "AGE", "SEX", "LINE", "COMP", input.gene()]).to_pandas()
        data = pd.merge(rpkm_corr, gene_corr, on = ['gene', 'AGE', 'SEX', 'LINE'])
        data["AGE"] = pd.Categorical(data["AGE"], categories=ages, ordered=True)
        data["SEX"] = pd.Categorical(data["SEX"], categories=sexes, ordered=True)
        data["LINE"] = pd.Categorical(data["LINE"], categories=lines, ordered=True)
        data["COMP"] = pd.Categorical(data["COMP"], categories=comp_order, ordered=True)
        sorted_data = data.sort_values(['LINE', 'AGE', 'SEX', 'COMP'])
        match input.filter():
            case "1":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData["GROUP"] = outputData["AGE"]
                return outputData
            case "2":
                outputData = sorted_data.loc[sorted_data["SEX"].isin(input.sex()), ("SEX", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData["GROUP"] = outputData["SEX"]
                return outputData
            case "3":
                outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData["GROUP"] = outputData["LINE"]
                return outputData
            case "4":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
                return outputData
            case "5":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "6":
                outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("SEX", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "7":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
                return outputData

    @render_widget
    def gene_corr_plot():
        if input.gene() not in sorted_genes:
            return
        fig = make_subplots(rows=1, cols=3)
        corr_data = filtered_gene_corr()
        position = 1

        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"

        legendKey = [False, False, True]

        for comp in corr_data['COMP'].unique():
            comp_data = corr_data[corr_data['COMP'] == comp]

            # ---- Scatter points for each group ----
            for group in comp_data['GROUP'].unique():
                subset = comp_data[comp_data['GROUP'] == group]

                x = subset[f"{input.gene()}_x"]
                y = subset[f"{input.gene()}_y"]

                # Remove NaNs and ensure x>0 for log scale
                mask = x.notna() & y.notna() & (x > 0)
                x = x[mask]
                y = y[mask]

                # Floor y values at 0
                y = np.maximum(y, 0)

                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        name=group,
                        marker_color=color_map[group],
                        mode="markers",
                        legendgroup="test",
                        showlegend=legendKey[position-1]
                    ),
                    row=1,
                    col=position
                )

            # ---- One regression line across all groups ----
            x_all = comp_data[f"{input.gene()}_x"]
            y_all = comp_data[f"{input.gene()}_y"]

            mask_all = x_all.notna() & y_all.notna() & (x_all > 0)
            x_all = x_all[mask_all]
            y_all = np.maximum(y_all[mask_all], 0)

            if len(x_all) > 1:
                x_log = np.log10(x_all)
                result = linregress(x_log, y_all)

                x_sorted = np.sort(x_all)
                y_fit = result.slope * np.log10(x_sorted) + result.intercept
                y_fit = np.maximum(y_fit, 0)  # keep regression line >=0

                fig.add_trace(
                    go.Scatter(
                        x=x_sorted,
                        y=y_fit,
                        mode="lines",
                        line=dict(color="black", width=2),
                        name="Trendline",
                        showlegend=False
                    ),
                    row=1,
                    col=position
                )

            position += 1

        fig.update_layout(
            title="Gene Body Correlation: " + input.gene(),
            showlegend=True,
            title_font=dict(size=24, weight="bold", color=fontcolor),
            font=dict(color=fontcolor),
            paper_bgcolor=papercolor,
            plot_bgcolor=bgcolor,
            margin=dict(t=100)
        )

        fig.update_yaxes(title_text="5modCG (%)", row=1, col=1)
        fig.update_yaxes(title_text="5mCG (%)", row=1, col=2)
        fig.update_yaxes(title_text="5hmCG (%)", row=1, col=3)

        fig.update_xaxes(type="log", title_text="log (RPKM)", row=1, col=1)
        fig.update_xaxes(type="log", title_text="log (RPKM)", row=1, col=2)
        fig.update_xaxes(type="log", title_text="log (RPKM)", row=1, col=3)

        return fig
    
    @reactive.Calc
    def filtered_tss_corr() -> pd.DataFrame:
        rpkm_corr = pq.read_table('ALL_RPKM_DATA_FILTERED_T_v2.parquet', columns=["gene", "AGE", "SEX", "LINE", input.gene()]).to_pandas()
        tss_corr = pq.read_table('ALL_TSS_PER_SAMPLE_T_v2.parquet', columns=["gene", "AGE", "SEX", "LINE", "COMP", input.gene()]).to_pandas()
        data = pd.merge(rpkm_corr, tss_corr, on = ['gene', 'AGE', 'SEX', 'LINE'])
        data["AGE"] = pd.Categorical(data["AGE"], categories=ages, ordered=True)
        data["SEX"] = pd.Categorical(data["SEX"], categories=sexes, ordered=True)
        data["LINE"] = pd.Categorical(data["LINE"], categories=lines, ordered=True)
        data["COMP"] = pd.Categorical(data["COMP"], categories=comp_order, ordered=True)
        sorted_data = data.sort_values(['LINE', 'AGE', 'SEX', 'COMP'])
        match input.filter():
            case "1":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData["GROUP"] = outputData["AGE"]
                return outputData
            case "2":
                outputData = sorted_data.loc[sorted_data["SEX"].isin(input.sex()), ("SEX", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData["GROUP"] = outputData["SEX"]
                return outputData
            case "3":
                outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData["GROUP"] = outputData["LINE"]
                return outputData
            case "4":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
                return outputData
            case "5":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "6":
                outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("SEX", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData["GROUP"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
                return outputData
            case "7":
                outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", f"{input.gene()}_x", f"{input.gene()}_y", "COMP")]
                outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
                outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
                outputData["GROUP"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
                return outputData

    @render_widget
    def tss_corr_plot():
        if input.gene() not in sorted_genes:
            return
        fig = make_subplots(rows=1, cols=3)
        corr_data = filtered_tss_corr()
        position = 1

        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"

        legendKey = [False, False, True]
        
        for comp in corr_data['COMP'].unique():
            comp_data = corr_data[corr_data['COMP'] == comp]

            # ---- Scatter points for each group ----
            for group in comp_data['GROUP'].unique():
                subset = comp_data[comp_data['GROUP'] == group]

                x = subset[f"{input.gene()}_x"]
                y = subset[f"{input.gene()}_y"]

                # Remove NaNs and ensure x>0 for log scale
                mask = x.notna() & y.notna() & (x > 0)
                x = x[mask]
                y = y[mask]

                # Floor y values at 0
                y = np.maximum(y, 0)

                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        name=group,
                        marker_color=color_map[group],
                        mode="markers",
                        legendgroup="test",
                        showlegend=legendKey[position-1]
                    ),
                    row=1,
                    col=position
                )

            # ---- One regression line across all groups ----
            x_all = comp_data[f"{input.gene()}_x"]
            y_all = comp_data[f"{input.gene()}_y"]

            mask_all = x_all.notna() & y_all.notna() & (x_all > 0)
            x_all = x_all[mask_all]
            y_all = np.maximum(y_all[mask_all], 0)

            if len(x_all) > 1:
                x_log = np.log10(x_all)
                result = linregress(x_log, y_all)

                x_sorted = np.sort(x_all)
                y_fit = result.slope * np.log10(x_sorted) + result.intercept
                y_fit = np.maximum(y_fit, 0)  # keep regression line >=0

                fig.add_trace(
                    go.Scatter(
                        x=x_sorted,
                        y=y_fit,
                        mode="lines",
                        line=dict(color="black", width=2),
                        name="Trendline",
                        showlegend=False
                    ),
                    row=1,
                    col=position
                )

            position += 1

        fig.update_layout(
            title="Promoter Correlation: " + input.gene(),
            showlegend=True,
            title_font=dict(size=24, weight="bold", color=fontcolor),
            font=dict(color=fontcolor),
            paper_bgcolor=papercolor,
            plot_bgcolor=bgcolor,
            margin=dict(t=100)
        )

        fig.update_yaxes(title_text="5modCG (%)", row=1, col=1)
        fig.update_yaxes(title_text="5mCG (%)", row=1, col=2)
        fig.update_yaxes(title_text="5hmCG (%)", row=1, col=3)

        fig.update_xaxes(type="log", title_text="log (RPKM)", row=1, col=1)
        fig.update_xaxes(type="log", title_text="log (RPKM)", row=1, col=2)
        fig.update_xaxes(type="log", title_text="log (RPKM)", row=1, col=3)

        return fig
    
app = App(app_ui, server)