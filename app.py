import pyarrow.parquet as pq
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shinywidgets import output_widget, render_widget

genes = pd.read_csv("DETECTED_GENES.csv", header=None, index_col=False)

ages = ["Young", "Adult", "Old"]
sexes = ["Female", "Male"]
lines = ["Astrocytes", "Neurons", "Microglia"]
comp_order = ["modCG", "mCG", "hmCG"]
sorted_genes = sorted(genes[0])
color_map = {
    'Young': "#74a5ce",
    'Adult': "#afc75b",
    'Old': "#b14646",
    'Male': "#4c3ec7",
    'Female': "#b33ba9",
    'Astrocytes': "#4d26bb",
    'Neurons': "#53410F",
    'Microglia': "#bb58da",
    'Young Male': "#4ecbeb",
    'Young Female': "#f088e7",
    'Adult Male': "#6d93e6",
    'Adult Female': "#c65cdb",
    'Old Male': "#4224ee",
    'Old Female': "#eb25eb",
    'Young Astrocytes': "#9694f3",
    'Young Neurons': "#ecd279",
    'Young Microglia': "#d894f3",
    'Adult Astrocytes': "#6b58d4",
    'Adult Neurons': "#ebb238",
    'Adult Microglia': "#bc58d4",
    'Old Astrocytes': "#130ff3",
    'Old Neurons': "#b37708",
    'Old Microglia': "#cd0ff3",
    'Male Astrocytes': "#0e49eb",
    'Female Astrocytes': "#7708f5",
    'Male Neurons': "#e45311",
    'Female Neurons': "#e9a617",
    'Male Microglia': "#a351f0",
    'Female Microglia': "#dc24ec",
    'Young Female Astrocytes': "#eca57c",
    'Young Female Neurons': "#f0d68e",
    'Young Female Microglia': "#e5baee",
    'Adult Female Astrocytes': "#eb8d40",
    'Adult Female Neurons': "#926C3B",
    'Adult Female Microglia': "#a550c7",
    'Old Female Astrocytes': "#e44040",
    'Old Female Neurons': "#5f3404",
    'Old Female Microglia': "#4c0d5f",
    'Young Male Astrocytes': "#97b8e2",
    'Young Male Neurons': "#8ad6c6",
    'Young Male Microglia': "#89d8a7",
    'Adult Male Astrocytes': "#49b3e4",
    'Adult Male Neurons': "#319988",
    'Adult Male Microglia': "#57A063",
    'Old Male Astrocytes': "#4839c9",
    'Old Male Neurons': "#094250",
    'Old Male Microglia': "#043608"
}

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_selectize(
            "gene", 
            "Gene", 
            sorted_genes,
            selected="Shh"
        ),
        ui.input_select(
            "filter", 
            "Grouping", 
            {
                1: "Age",
                2: "Sex",
                3: "Cell Line",
                4: "Age and Sex",
                5: "Age and Cell Line",
                6: "Sex and Cell Line",
                7: "Age, Sex, and Cell Line"
            },
            selected=7
        ),
        ui.panel_conditional(
            "input.filter == 1 || input.filter == 4 || input.filter == 5 || input.filter == 7",
            ui.input_checkbox_group(
                "age",
                "Filter by Age",
                ["Young", "Adult", "Old"],
                selected=ages
            ),
        ),
        ui.panel_conditional(
            "input.filter == 2 || input.filter == 4 || input.filter == 6 || input.filter == 7",
            ui.input_checkbox_group(
                "sex",
                "Filter by Sex",
                ["Female", "Male"],
                selected=sexes
            ),
        ),
        ui.panel_conditional(
            "input.filter == 3 || input.filter == 5 || input.filter == 6 || input.filter == 7",
            ui.input_checkbox_group(
                "line",
                "Filter by Cell Line",
                ["Astrocytes", "Neurons", "Microglia"],
                selected=lines
            ),
        ),
        ui.download_button("download_data", "Download Data"),
        ui.download_button("download_plot", "Download Plot")
    ),
    ui.page_navbar(
        ui.nav_panel(
            "Gene Expression",
            ui.layout_columns(
                output_widget("expression_plot")
            )
        ),
        ui.nav_panel(
            "Gene Methylation",
            ui.layout_columns(
                output_widget("gene_body_plot")
            ),
            ui.layout_columns(
                output_widget("tss_plot")
            )
        )
    )
)

def server(input: Inputs, output: Outputs, session: Session):
    @render.download(filename="out_data.csv")
    def download_data():
        outData = filtered_expr()
        return outData.to_csv("out_data.csv", index=False)
    
    @render.download(filename="out_plot.png")
    def download_plot():
        outPlot = expression_plot()
        return outPlot.write_image("out_plot.png", index=False)
    
    @reactive.Calc
    def filtered_expr() -> pd.DataFrame:
        data = pq.read_table('ALL_RPKM_LABELED_FILTERED.parquet', columns=["AGE", "SEX", "LINE", input.gene()]).to_pandas()
        data["AGE"] = pd.Categorical(data["AGE"], categories=ages, ordered=True)
        data["SEX"] = pd.Categorical(data["SEX"], categories=sexes, ordered=True)
        data["LINE"] = pd.Categorical(data["LINE"], categories=lines, ordered=True)
        sorted_data = data.sort_values(['LINE', 'SEX', 'AGE'])
        if input.filter() == "1":
            outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "LINE", input.gene())]
            return outputData
        if input.filter() == "2":
            outputData = sorted_data.loc[sorted_data["SEX"].isin(input.sex()), ("SEX", "LINE", input.gene())]
            return outputData
        if input.filter() == "3":
            outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("LINE", input.gene())]
            return outputData
        if input.filter() == "4":
            outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData["AGE_SEX"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
            return outputData
        if input.filter() == "5":
            outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "LINE", input.gene())]
            outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
            outputData["AGE_LINE"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
            return outputData
        if input.filter() == "6":
            outputData = sorted_data.loc[sorted_data["LINE"].isin(input.line()), ("SEX", "LINE", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData["LINE_SEX"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
            return outputData
        if input.filter() == "7":
            outputData = sorted_data.loc[sorted_data["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
            outputData["ALL"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
            return outputData

    @render_widget
    def expression_plot():
        if input.gene() not in genes:
            return
        fig = go.Figure()
        data = filtered_expr()
        if input.filter() == "1":
            for group in data['AGE'].unique():
                fig.add_trace(go.Box(y = data[data['AGE'] == group][input.gene()],
                                     boxpoints = 'all', jitter = 0.5, 
                                     pointpos = 0, name = group, marker_color=color_map[group]))
        elif input.filter() == "2":
            for group in data['SEX'].unique():
                fig.add_trace(go.Box(y = data[data['SEX'] == group][input.gene()],
                                     boxpoints = 'all', jitter = 0.5, 
                                     pointpos = 0, name = group, marker_color=color_map[group]))
        elif input.filter() == "3":
            for group in data['LINE'].unique():
                fig.add_trace(go.Box(y = data[data['LINE'] == group][input.gene()],
                                     boxpoints = 'all', jitter = 0.5, 
                                     pointpos = 0, name = group, marker_color=color_map[group]))
        elif input.filter() == "4":
            for group in data['AGE_SEX'].unique():
                fig.add_trace(go.Box(y = data[data['AGE_SEX'] == group][input.gene()],
                                     boxpoints = 'all', jitter = 0.5, 
                                     pointpos = 0, name = group, marker_color=color_map[group]))
        elif input.filter() == "5":
            for group in data['AGE_LINE'].unique():
                fig.add_trace(go.Box(y = data[data['AGE_LINE'] == group][input.gene()],
                                     boxpoints = 'all', jitter = 0.5, 
                                     pointpos = 0, name = group, marker_color=color_map[group]))
        elif input.filter() == "6":
            for group in data['LINE_SEX'].unique():
                fig.add_trace(go.Box(y = data[data['LINE_SEX'] == group][input.gene()],
                                     boxpoints = 'all', jitter = 0.5, 
                                     pointpos = 0, name = group, marker_color=color_map[group]))
        elif input.filter() == "7":
            for group in data['ALL'].unique():
                y = data[data['ALL'] == group][input.gene()]
                fig.add_trace(go.Box(y = data[data['ALL'] == group][input.gene()],
                                     boxpoints = 'all', jitter = 0.5, 
                                     pointpos = 0, name = group, marker_color=color_map[group]))
        fig.update_layout(
            title = input.gene(),
            title_font = dict(
                size = 24,
                textcase = "upper",
                weight = "bold"
            ),
            yaxis_title = "RPKM",
            xaxis_title = "",
            showlegend = False,
            plot_bgcolor = "#f0f0f0"
        )
        return fig
    
    @reactive.Calc
    def filtered_body() -> pd.DataFrame:
        gene_corr = pq.read_table('ALL_GENE_BODY_PER_SAMPLE_TRANSPOSED_FINAL.parquet', columns=["AGE", "SEX", "LINE", "COMP", input.gene()]).to_pandas()
        gene_corr["AGE"] = pd.Categorical(gene_corr["AGE"], categories=ages, ordered=True)
        gene_corr["SEX"] = pd.Categorical(gene_corr["SEX"], categories=sexes, ordered=True)
        gene_corr["LINE"] = pd.Categorical(gene_corr["LINE"], categories=lines, ordered=True)
        gene_corr["COMP"] = pd.Categorical(gene_corr["COMP"], categories=comp_order, ordered=True)
        sorted_body = gene_corr.sort_values(['LINE', 'SEX', 'AGE', 'COMP'])
        if input.filter() == "1":
            outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
            return outputData
        if input.filter() == "2":
            outputData = sorted_body.loc[sorted_body["SEX"].isin(input.sex()), ("SEX", "LINE", "COMP", input.gene())]
            return outputData
        if input.filter() == "3":
            outputData = sorted_body.loc[sorted_body["LINE"].isin(input.line()), ("LINE", "COMP", input.gene())]
            return outputData
        if input.filter() == "4":
            outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData["AGE_SEX"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
            return outputData
        if input.filter() == "5":
            outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
            outputData["AGE_LINE"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
            return outputData
        if input.filter() == "6":
            outputData = sorted_body.loc[sorted_body["LINE"].isin(input.line()), ("SEX", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData["LINE_SEX"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
            return outputData
        if input.filter() == "7":
            outputData = sorted_body.loc[sorted_body["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
            outputData["ALL"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
            return outputData

    @render_widget
    def gene_body_plot():
        if input.gene() not in genes:
            return
        fig = make_subplots(rows=1, cols=3)
        body_data = filtered_body()
        position=1
        for comp in body_data['COMP'].unique():
            comp_data = body_data[body_data['COMP'] == comp]
            if input.filter() == "1":
                for group in comp_data['AGE'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['AGE'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "2":
                for group in comp_data['SEX'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['SEX'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "3":
                for group in comp_data['LINE'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['LINE'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "4":
                for group in comp_data['AGE_SEX'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['AGE_SEX'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "5":
                for group in comp_data['AGE_LINE'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['AGE_LINE'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "6":
                for group in comp_data['LINE_SEX'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['LINE_SEX'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "7":
                for group in comp_data['ALL'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['ALL'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            position += 1
        fig.update_layout(
            title = "Gene Body Methylation",
            showlegend = False,
            plot_bgcolor = "#f0f0f0"
        )

        fig.update_yaxes(title_text="modCG (%)", range=[0, 100], row=1, col=1)
        fig.update_yaxes(title_text="mCG (%)", range=[0, 100], row=1, col=2)
        fig.update_yaxes(title_text="hmCG (%)", range=[0, 100], row=1, col=3)
        return fig

    @reactive.Calc
    def filtered_tss() -> pd.DataFrame:
        tss_corr = pq.read_table('ALL_TSS_PER_SAMPLE_TRANSPOSED.parquet', columns=["AGE", "SEX", "LINE", "COMP", input.gene()]).to_pandas()
        tss_corr["AGE"] = pd.Categorical(tss_corr["AGE"], categories=ages, ordered=True)
        tss_corr["SEX"] = pd.Categorical(tss_corr["SEX"], categories=sexes, ordered=True)
        tss_corr["LINE"] = pd.Categorical(tss_corr["LINE"], categories=lines, ordered=True)
        tss_corr["COMP"] = pd.Categorical(tss_corr["COMP"], categories=comp_order, ordered=True)
        sorted_tss = tss_corr.sort_values(['LINE', 'SEX', 'AGE', 'COMP'])
        if input.filter() == "1":
            outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
            return outputData
        if input.filter() == "2":
            outputData = sorted_tss.loc[sorted_tss["SEX"].isin(input.sex()), ("SEX", "LINE", "COMP", input.gene())]
            return outputData
        if input.filter() == "3":
            outputData = sorted_tss.loc[sorted_tss["LINE"].isin(input.line()), ("LINE", "COMP", input.gene())]
            return outputData
        if input.filter() == "4":
            outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData["AGE_SEX"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)
            return outputData
        if input.filter() == "5":
            outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
            outputData["AGE_LINE"] = outputData["AGE"].astype(str) + " " + outputData["LINE"].astype(str)
            return outputData
        if input.filter() == "6":
            outputData = sorted_tss.loc[sorted_tss["LINE"].isin(input.line()), ("SEX", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData["LINE_SEX"] = outputData["SEX"].astype(str) + " " + outputData["LINE"].astype(str)
            return outputData
        if input.filter() == "7":
            outputData = sorted_tss.loc[sorted_tss["AGE"].isin(input.age()), ("AGE", "SEX", "LINE", "COMP", input.gene())]
            outputData = outputData.loc[outputData["SEX"].isin(input.sex()), ]
            outputData = outputData.loc[outputData["LINE"].isin(input.line()), ]
            outputData["ALL"] = outputData["AGE"].astype(str) + " " + outputData["SEX"].astype(str)  + " " + outputData["LINE"].astype(str) 
            return outputData

    @render_widget
    def tss_plot():
        if input.gene() not in genes:
            return
        fig = make_subplots(rows=1, cols=3)
        tss_data = filtered_tss()
        position=1
        for comp in tss_data['COMP'].unique():
            comp_data = tss_data[tss_data['COMP'] == comp]
            if input.filter() == "1":
                for group in comp_data['AGE'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['AGE'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "2":
                for group in comp_data['SEX'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['SEX'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "3":
                for group in comp_data['LINE'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['LINE'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "4":
                for group in comp_data['AGE_SEX'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['AGE_SEX'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "5":
                for group in comp_data['AGE_LINE'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['AGE_LINE'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "6":
                for group in comp_data['LINE_SEX'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['LINE_SEX'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            elif input.filter() == "7":
                for group in comp_data['ALL'].unique():
                    fig.add_trace(go.Box(y = comp_data[comp_data['ALL'] == group][input.gene()],
                                        boxpoints = 'all', jitter = 0.5, 
                                        pointpos = 0, name = group, marker_color=color_map[group]),
                                row=1, col=position)
            position += 1
        fig.update_layout(
            title = "Promoter Methylation",
            showlegend = False,
            plot_bgcolor = "#f0f0f0"
        )
        
        fig.update_yaxes(title_text="modCG (%)", range=[0, 100], row=1, col=1)
        fig.update_yaxes(title_text="mCG (%)", range=[0, 100], row=1, col=2)
        fig.update_yaxes(title_text="hmCG (%)", range=[0, 100], row=1, col=3)
        return fig

app = App(app_ui, server)