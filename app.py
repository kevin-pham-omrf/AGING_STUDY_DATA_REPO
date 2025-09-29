import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shinywidgets import output_widget, render_widget

data = pd.read_csv("ALL_RPKM_LABELED_SUBSET.csv", na_values = "NA", header=0)

genes = data.select_dtypes(include=["float64"]).columns.tolist()
ages = data["AGE"].unique().tolist()
sexes = data["SEX"].unique().tolist()
lines = data["LINE"].unique().tolist()

age_order = ["Young", "Adult", "Old"]
sex_order = ["Female", "Male"]
line_order = ["Astrocytes", "Neurons", "Microglia"]

data["AGE"] = pd.Categorical(data["AGE"], categories=age_order, ordered=True)
data["SEX"] = pd.Categorical(data["SEX"], categories=sex_order, ordered=True)
data["LINE"] = pd.Categorical(data["LINE"], categories=line_order, ordered=True)

sorted_data = data.sort_values(['LINE', 'SEX', 'AGE'])

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_selectize(
            "gene", 
            "Gene", 
            sorted(genes),
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
    output_widget("boxplot")
)

def server(input: Inputs, output: Outputs, session: Session):
    @render.download(filename="out_data.csv")
    def download_data():
        outData = filtered_df()
        return outData.to_csv("out_data.csv", index=False)
    
    @render.download(filename="out_plot.png")
    def download_plot():
        outPlot = boxplot()
        return outPlot.write_image("out_plot.png", index=False)
    
    @reactive.Calc
    def filtered_df() -> pd.DataFrame:
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
    def boxplot():
        fig = go.Figure()
        data = filtered_df()
        if input.gene() not in genes:
            return
        if input.filter() == "1":
            for group in data['AGE'].unique():
                fig.add_trace(go.Box(y = data[data['AGE'] == group][input.gene()],
                                     boxpoints='all', jitter=0.5, pointpos=0, name=group))
        elif input.filter() == "2":
            for group in data['SEX'].unique():
                fig.add_trace(go.Box(y = data[data['SEX'] == group][input.gene()],
                                     boxpoints='all', jitter=0.5, pointpos=0, name=group))
        elif input.filter() == "3":
            for group in data['LINE'].unique():
                fig.add_trace(go.Box(y = data[data['LINE'] == group][input.gene()],
                                     boxpoints='all', jitter=0.5, pointpos=0, name=group))
        elif input.filter() == "4":
            for group in data['AGE_SEX'].unique():
                fig.add_trace(go.Box(y = data[data['AGE_SEX'] == group][input.gene()],
                                     boxpoints='all', jitter=0.5, pointpos=0, name=group))
        elif input.filter() == "5":
            for group in data['AGE_LINE'].unique():
                fig.add_trace(go.Box(y = data[data['AGE_LINE'] == group][input.gene()],
                                     boxpoints='all', jitter=0.5, pointpos=0, name=group))
        elif input.filter() == "6":
            for group in data['LINE_SEX'].unique():
                fig.add_trace(go.Box(y = data[data['LINE_SEX'] == group][input.gene()],
                                     boxpoints='all', jitter=0.5, pointpos=0, name=group))
        elif input.filter() == "7":
            for group in data['ALL'].unique():
                fig.add_trace(go.Box(y = data[data['ALL'] == group][input.gene()],
                                     boxpoints='all', jitter=0.5, pointpos=0, name=group))
        fig.update_layout(
            title=input.gene(),
            title_font=dict(
                size=24,
                textcase="upper",
                weight="bold"
            ),
            yaxis_title="RPKM",
            xaxis_title="",
            showlegend=False,
            plot_bgcolor="#f0f0f0"
        )
        return fig

app = App(app_ui, server)