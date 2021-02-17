from dash import Dash, no_update
import dash_html_components as html
import dash_core_components as dcc
import dash_table, dash_bio
from dash.dependencies import Input, Output, State

import os
import pandas as pd

import dash_cytoscape as cyto
cyto.load_extra_layouts()

app = Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server


# Define Helper Functions

def table_data(df):
    return df[['name', 'length'] + [name for name in df.columns if name.endswith("score")]]

# Load Data
dataset_paths = []
dir_path = os.path.dirname(os.path.realpath(__file__))
df = pd.read_csv("/home/abc/Workspace/ConservedUTRs/alignments/three_alignment_table (copy).csv")
df["average_score"] = df[[name for name in df.columns if name.endswith("score")]].mean(axis=1)
df = df.sort_values(by=["average_score"], ascending=False).round(decimals=3)
data = ">No_Data\n"

#Layout
app.layout = html.Div(
    id="root",
    children=[
        html.Div(
            id="app-container",
            children=[
                
                html.Div(
                    id='tabs-container',
                    className="container",
                    children=[
                        html.Div(
                            id="headline",
                            children=[
                                html.H3(children="Conserved UTRs"),
                            ],
                        ),
                        dcc.Tabs(
                            id='data-tabs', 
                            value='table', 
                            children=[
                                dcc.Tab(
                                    id="table-tab",
                                    className='custom-tab',
                                    selected_className='custom-tab--selected',
                                    label='Table',
                                    value='table',
                                    children=[
                                        html.Div(
                                            id="table-container",
                                            className="container",
                                            children=[
                                                dash_table.DataTable(
                                                    id='table',
                                                    columns=[{"name": i[:-6], "id": i} if  i.endswith("score") else {"name": i, "id": i} for i in table_data(df).columns],
                                                    data=table_data(df).to_dict('records'),
                                                    sort_action="native"
                                                )
                                            ]
                                        ),
                                        html.Div(
                                            id="alignment-container",
                                            children=[
                                                dash_bio.AlignmentChart(
                                                    id='my-alignment-viewer',
                                                    data=data,
                                                    showconservation=False,
                                                    showgap=False,
                                                    colorscale={'A':'#1111ff', 'T':'#1111aa', 'C':'#ff1111', 'G':'#aa1111', '-': '#ffffff'},
                                                    overview='none',
                                                    showconsensus=False,
                                                    height=110
                                                )
                                            ]
                                        ),
                                        html.Div(
                                            className="controls-block",
                                            children=[
                                                html.Div(
                                                    className="control-element",
                                                    children=[
                                                        dcc.Checklist(
                                                            id='species-checklist',
                                                            options=[
                                                                {'label': i[:-6], 'value': i} for i in table_data(df).columns if (i.endswith("score") and not i.startswith("average"))
                                                            ],
                                                            inputClassName = "padded"
                                                        )
                                                    ]
                                                )
                                            ]
                                        ),
                                    ]
                                ),
                                dcc.Tab(
                                    id="about-tab",
                                    className='custom-tab',
                                    selected_className='custom-tab--selected',
                                    label='About',
                                    value='about',
                                    children=[
                                        "teesetset"
                                    ]
                                )
                            ]
                        )
                    ]
                )           
            ]
        )
    ]
)

@app.callback(Output('my-alignment-viewer', 'data'),
              Input('table', 'active_cell'))
def get_active_alignment(active_cell):
    if active_cell:
        if active_cell['column_id'] not in ['name', 'length', 'average_score']:
            my_fasta = ""
            my_fasta += ">{}\n{}-\n".format(df.iloc[active_cell['row']]['name'], df.iloc[active_cell['row']][active_cell['column_id'][0:-6] + "_seqaln"])
            my_fasta += ">{}\n{}-".format(active_cell['column_id'][0:-6], df.iloc[active_cell['row']][active_cell['column_id'][0:-6] + "_refaln"])
            return my_fasta
        elif active_cell['column_id'] in ['name', 'average_score']:
            pass
    return no_update
    

if __name__ == '__main__':
    app.run_server(debug=True,port=8083,host='0.0.0.0');
