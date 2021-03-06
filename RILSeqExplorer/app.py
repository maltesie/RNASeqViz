from dash import Dash, no_update
import dash_html_components as html
import dash_core_components as dcc
import dash_table
import dash_bio as dashbio
from dash.dependencies import Input, Output, State

import os
import pandas as pd
import networkx as nx
import json
import numpy as np
import webbrowser

import dash_cytoscape as cyto
cyto.load_extra_layouts()

app = Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server

circos_graph_data = [
    {
      "color": "#ff5722",
      "source": {
        "id": "chr1",
        "start": 1186054,
        "end": 1478117
      },
      "target": {
        "id": "chr1",
        "start": 1550000,
        "end": 1578117
      }
    },
    {
      "color": "#ff5722",
      "source": {
        "id": "chr1",
        "start": 1000000,
        "end": 1100000
      },
      "target": {
        "id": "chr2",
        "start": 100000,
        "end": 200000
      }
    }]

# Define Helper Functions
def filter_threshold(df, threshold):
    return df[np.log(df['nb_fragments']) > threshold]

def filter_targets(df, node):
    return df[df['gene1'].str.match(node) | df['gene2'].str.match(node)]

def filter_cc(df, node):
    g = nx.Graph()
    g.add_edges_from(zip(df['gene1'], df['gene2']))
    for component in nx.connected_components(g):
        if node in component: break
    search_str = r"{}".format('|'.join(component))
    return df[df['gene1'].str.match(search_str) | df['gene2'].str.match(search_str)]

def filter_search(df, search_string):
    return df[df['name1'].str.contains(search_string, case=False) | df['name2'].str.contains(search_string, case=False)]

def table_data(df):
    return df[['name1', 'type1', 'name2', 'type2', 'libs', 'nb_fragments']]

def fix_ig_label(label):
    if '|' in label:
        return '\n'.join(['IG']+label[3:].split('|'))
    else: return label
    
def translate(value, leftMin, leftMax, rightMin, rightMax):
    leftSpan = leftMax - leftMin
    rightSpan = rightMax - rightMin
    valueScaled = min(max(float(value - leftMin) / float(leftSpan), 0.0), 1.0)
    return rightMin + (valueScaled * rightSpan)

def cb_edge_percentage(value):
    return int(round(translate(value, 0, 1, 18, 76)))

def cb_node_percentage(value):
    return int(round(translate(value, 0, 0.1, 14, 72)))

def cytoscape_data(df, norm):
    jsondata = {"nodes": [], "edges": []}
    nodes = jsondata["nodes"]
    edges = jsondata["edges"]
    fragment_count = {}
    current_norm = df['nb_fragments'].sum()
    types = {}
    
    for i, interaction in df.iterrows():
        source, target, type_source, type_target, relpos1, relpos2, nb_fragments, from_region, to_region = \
            interaction[['gene1', 'gene2', 'type1', 'type2', 'relpos1', 'relpos2', 'nb_fragments', 'name1', 'name2']]
        nb_fragments /= norm
        pos = max(relpos1, relpos2)
        if ((type_target == "srna") != (type_source=="srna")): edge_type = "srna_edge"
        else: edge_type = "other_edge"
        #source = fix_ig_label(source)
        #target = fix_ig_label(target)
        edges.append({'data': dict(source=source, target=target, fragments=nb_fragments, 
                                   from_region=from_region, to_region=to_region, pos=pos, typ=edge_type),
                      'classes': edge_type}) 
        if source in fragment_count: 
            fragment_count[source] += nb_fragments
        else: 
            fragment_count[source] = nb_fragments
            types[source] = type_source
        
        if target in fragment_count: 
            fragment_count[target] += nb_fragments
        else: 
            fragment_count[target] = nb_fragments
            types[target] = type_target
            
    for node in fragment_count:
        nodes.append({'data':dict(id=node, name=fix_ig_label(node), fragments=fragment_count[node], current_fragments=fragment_count[node]*norm/current_norm, typ=types[node]), 'classes':types[node]})
    
    return jsondata

# Load Data
dataset_paths = []
dir_path = os.path.dirname(os.path.realpath(__file__))
for file in os.listdir(os.path.join(dir_path, "assets")):
    if file.endswith(".csv"):
        dataset_paths.append(os.path.abspath(os.path.join(dir_path, "assets", file)))
        
with open(os.path.join(dir_path, "assets", 'mystylesheet.json')) as json_file:
    stylesheet = json.load(json_file)
    
initial_df = pd.read_csv(dataset_paths[0])
filtered_df = pd.read_csv(dataset_paths[0])
selected_node = None
min_edge, max_edge = np.log(initial_df['nb_fragments'].min()), np.log(initial_df['nb_fragments'].max())
fragments_sum = initial_df['nb_fragments'].sum()
elements = cytoscape_data(initial_df[:100], fragments_sum)

app.layout = html.Div(
    id="root",
    children=[
        html.Div(
            id="app-container",
            children=[
                html.Div(
                    id="left-column",
                    children=[
                        html.Div(
                            id="headline",
                            children=[
                                html.H3(children="RIL-Seq Analysis Visualization"),
                            ],
                        ),
                        html.Div(
                            className="container",
                            children=[
                                html.Div(
                                    className="controls-block horizontal",
                                    children=[
                                        html.Div(
                                            className="control-element",
                                            children=[
                                                html.P("Dataset:"),
                                                dcc.Dropdown(
                                                    id='dropdown-update-dataset',
                                                    value=dataset_paths[0],
                                                    clearable=False,
                                                    options=[
                                                        {'label': name.split('/')[-1], 'value': name}
                                                        for name in dataset_paths
                                                    ]
                                                ),
                                            ]
                                        ),
                                        html.Div(
                                            className="control-element",
                                            children=[
                                                html.P("Positioning:"),
                                                dcc.Dropdown(
                                                    id='dropdown-update-layout',
                                                    value='random',
                                                    clearable=False,
                                                    options=[
                                                        {'label': name.capitalize(), 'value': name}
                                                        for name in ['cose-bilkent', 'grid', 'random', 'circle', 'concentric']
                                                    ]
                                                )
                                            ]
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="controls-block",
                                    children=[
                                        html.Div(
                                            children=[
                                                #html.P("min reads:", className="line-label"),
                                                dcc.Slider(
                                                    id='reads-slider',
                                                    min=min_edge - 0.01*min_edge,
                                                    max=max_edge + 0.01*max_edge,
                                                    step=(max_edge-min_edge)/100,
                                                    value= min_edge + 0.5*(max_edge-min_edge),
                                                    marks={
                                                        min_edge: '{}'.format(int(np.exp(min_edge))),
                                                        min_edge + 0.5*(max_edge-min_edge): '{}'.format(int(np.exp((min_edge + 0.5*(max_edge-min_edge))))),
                                                        max_edge: '{}'.format(int(np.exp(max_edge)))
                                                    }
                                                ),
                                            ]
                                        ),
                                        dcc.RadioItems(
                                            id='filter-radio',
                                            options=[
                                                {'label': 'all interactions', 'value': 'all'},
                                                {'label': 'only targets of selected', 'value': 'targets'},
                                                {'label': 'connected component of selected', 'value': 'cc'}
                                            ],
                                            value='all'
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="controls-block horizontal",
                                    children=[
                                        html.P("Search:", className="line-label"),
                                        dcc.Input(
                                            id="input-search",
                                            type="text",
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="controls-block",
                                    children=[
                                        html.P(children=["No interaction selected."], id="test-output")
                                    ]
                                ),
                            ]
                        )
                    ],
                ),
                html.Div(
                    id='tabs-container',
                    className="container",
                    children=[
                        dcc.Tabs(
                            id='data-tabs', 
                            value='graph', 
                            children=[
                                dcc.Tab(
                                    id="graph-tab",
                                    className='custom-tab',
                                    selected_className='custom-tab--selected',
                                    label='Graph',
                                    value='graph',
                                    children=[
                                        html.Div(
                                            id="graph-container",
                                            className="container",
                                            children=[
                                                cyto.Cytoscape(
                                                    id='graph',
                                                    elements=elements,
                                                    stylesheet=stylesheet,
                                                    responsive=True,
                                                    layout={'name':'random'}
                                                ),
                                                html.Button('DOWNLOAD SVG', id='save-svg', n_clicks=0)
                                            ]
                                        ),
                                    ]
                                ),
                                dcc.Tab(
                                    id="circos-tab",
                                    className='custom-tab',
                                    selected_className='custom-tab--selected',
                                    label='Circos',
                                    value='circos',
                                    children=[
                                        html.Div(
                                            id="circos-container",
                                            className="container",
                                            children=[
                                                #html.Div(id='spacer1'),
                                                dashbio.Circos(
                                                    enableZoomPan=True,
                                                    enableDownloadSVG=True,
                                                    id='my-dashbio-circos',
                                                    config = {
                                                        'ticks': {
                                                            'display': True, 
                                                            'spacing': 100000,
                                                            'color': '#000',
                                                            'labelDenominator': 1000000,
                                                            'labelSuffix': ' Mb',
                                                            'majorSpacing':5,
                                                            'minorSpacing':1,
                                                            'labelSpacing':5
                                                            }
                                                        },
                                                    layout=[
                                                        {
                                                          "id": "chr1",
                                                          "label": "chr1",
                                                          "color": "#999999",
                                                          "len": 2961149
                                                        },
                                                        {
                                                          "id": "chr2",
                                                          "label": "chr2",
                                                          "color": "#CCCCCC",
                                                          "len": 1072315
                                                        }],
                                                    selectEvent={"0": "hover", "1": "click", "2": "both"},
                                                    tracks=[{
                                                        'type': 'CHORDS',
                                                        'data': circos_graph_data,
                                                        'config': {
                                                            'tooltipContent': {
                                                                'source': 'source',
                                                                'sourceID': 'id',
                                                                'target': 'target',
                                                                'targetID': 'id',
                                                                'targetEnd': 'end'
                                                            }
                                                        }
                                                    }]
                                                ),
                                                #html.Button('Save', id='save-circos', n_clicks=0)
                                            ]
                                        )
                                    ]
                                ),
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
                                                    columns=[{"name": i, "id": i} for i in table_data(initial_df).columns],
                                                    data=table_data(initial_df).to_dict('records')
                                                )
                                            ]
                                        )
                                    ]
                                ),
                                dcc.Tab(
                                    id="about-tab",
                                    className='custom-tab',
                                    selected_className='custom-tab--selected',
                                    label='About',
                                    value='about',
                                    children=[
                                        html.Div(
                                            id="about-container",
                                            className="container",
                                            children=[
                                                html.P("Color gradient explanation"),
                                                html.Div(
                                                    className="controls-block horizontal deflate",
                                                    children=[
                                                        html.P("5'UTR", className="cb-left"),
                                                        
                                                        html.Div(
                                                            id = "colorbar-edges",
                                                            className="controls-block",
                                                        ),
                                                        html.P("3'UTR", className="cb-right"),
                                                        html.P("Edges have two types. Edges of type 1 are between a sRNA and anything else. \
                                                               All other edges are of type 2. Edges of type 2 are colored grey and edges of \
                                                                   type 1 follow the color gradient on the left. The color is defined by \
                                                                       the average position of the interaction mapped on the length of the \
                                                                           gene that interacts with the sRNA.", className="text-block")
                                                    ]
                                                ),
                                                html.Div(
                                                    className="controls-block horizontal deflate",
                                                    children=[
                                                        html.P("0%", className="cb-left"),
                                                        
                                                        html.Div(
                                                            id = "colorbar-nodes",
                                                            className="controls-block",
                                                        ),
                                                        html.P("10%", className="cb-right"),
                                                        html.P("sRNA nodes are colored in turquoise. All other nodes follow a relative coloring \
                                                               scheme that is normalized over the number of reads in the selected, visible graph. \
                                                                   The scale starts at 0% and is capped at 10%.", className="text-block")
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
            ]
        )
    ]
)


@app.callback(
    [Output('table', 'data'),
    Output('graph', 'elements'),
    Output('reads-slider', 'marks')],
    [Input('reads-slider', 'value'),
    Input('filter-radio', 'value'),
    Input('input-search', 'value')])
def update_selected_data(slider_value, radio_filter, search_string):
    global filtered_df
    filtered_df = filter_threshold(initial_df, slider_value)
    if (radio_filter == 'targets') and (selected_node is not None): filtered_df = filter_targets(filtered_df, selected_node)
    elif (radio_filter == 'cc') and (selected_node is not None): filtered_df = filter_cc(filtered_df, selected_node)
    if (search_string is not None) and (len(search_string) > 2): filtered_df = filter_search(filtered_df, search_string)
    return table_data(filtered_df).to_dict('records'), cytoscape_data(filtered_df, fragments_sum), {slider_value: '{} reads'.format(int(np.exp(slider_value)))}

@app.callback(
    [Output('filter-radio', 'options'),
     Output('test-output', 'children'),
     Output('filter-radio', 'value')],
    [Input('graph', 'selectedNodeData'),
     Input('graph', 'selectedEdgeData')],
    state=State('filter-radio', 'value'))
def set_selected_element(node_data, edge_data, radio_value):
    global selected_node
    if (node_data is None) or (node_data == []):
        selected_node = None
        radio_return = [{'label': 'all interactions', 'value': 'all'},
                        {'label': 'only targets of selected RNA', 'value': 'targets'},
                        {'label': 'connected component of selected RNA', 'value': 'cc'}]
        text_return = ["Select a node or edge in the graph."]
        if radio_value != 'all': value_return = 'all'
        else: value_return = no_update
        if (edge_data is None) or (edge_data == []):
            text_return = ["Select an interaction to display further information."]
        else:  
            edge_data = edge_data[0]
            text_return = ["reads 1: {}".format(edge_data['source']), html.Br(), 
                           "reads 2: {}".format(edge_data['target']), html.Br(), 
                           "# of reads: {}".format(int(fragments_sum*float(edge_data['fragments'])))]
            value_return = no_update
            if edge_data['typ'] == 'srna_edge':
                text_return.append(html.Div(className="horizontal deflate", children=[html.P("5'UTR", className="cb-left"),
                                                                                        html.Div(
                                                                                            id = "colorbar-edges",
                                                                                            className="controls-block",
                                                                                            children=[
                                                                                                html.Div(
                                                                                                    className='colorbar-arrow',
                                                                                                    style={
                                                                                                        'left': '{}%'.format(cb_edge_percentage(edge_data['pos'])),
                                                                                                    }
                                                                                                )
                                                                                            ]
                                                                                        ),
                                                                                        html.P("3'UTR", className="cb-right")]))

    else: 
        node_data = node_data[0]
        selected_node = node_data["id"]
        radio_return = [{'label': 'all interactions', 'value': 'all'},
                        {'label': 'only targets of {}'.format(selected_node), 'value': 'targets'},
                        {'label': 'connected component of {}'.format(selected_node), 'value': 'cc'}]
        selected_node_interactions = filter_targets(filtered_df, selected_node)
        nb_targets = len(set(list(selected_node_interactions['gene1'])+list(selected_node_interactions['gene2']))) - 1
        text_return = ["{} is involved in {} interactions with {} targets.".format(selected_node, len(selected_node_interactions), nb_targets)]
        value_return = no_update
        if node_data['typ'] == 'gene':
            text_return.append(html.Div(className="horizontal deflate", children=[html.P("0%", className="cb-left"),
                                                                                    html.Div(
                                                                                        id = "colorbar-nodes",
                                                                                        className="controls-block",
                                                                                        children=[
                                                                                            html.Div(
                                                                                                className='colorbar-arrow',
                                                                                                style={
                                                                                                    'left': '{}%'.format(cb_node_percentage(selected_node_interactions['nb_fragments'].sum()/filtered_df['nb_fragments'].sum())),
                                                                                                }
                                                                                            )
                                                                                        ]
                                                                                    ),
                                                                                    html.P("10%", className="cb-right")]))
    return radio_return, text_return, value_return

@app.callback(
    Output('graph', 'layout'),
    Input('dropdown-update-layout', 'value'))
def update_layout(layout):
    if layout == 'cose-bilkent':
        return {
            'name':'cose-bilkent',
            'quality': 'draft',
            'idealEdgeLength': 50,
            'nodeOverlap': 0,
            'refresh': 10,
            'fit': True,
            'randomize': False,
            'componentSpacing': 5,
            'nodeRepulsion': 5000,
            'edgeElasticity': 0.5,
            'nestingFactor': 0.1,
            'gravity': 0.25,
            'numIter': 300,
            'gravityRange': 5,
            'animate': False
        }
    elif layout == 'concentric':
        return {
            'name':'concentric',
            'animate':False
        }
    else:
        return {
            'name': layout,
            'animate': False
        }
@app.callback(
    Output("graph", "generateImage"),
    Input('save-svg', 'n_clicks'))
def get_image(clicks):
    if clicks>0:
        return {
            'type': 'svg',
            ''
            'action': 'download'
            }
    else: return no_update
"""
@app.callback(
    Output("circos", "SVG"),
    Input('save-circos', 'n_clicks'))
def get_image_circos(clicks):
    if clicks>0:
        return {
            'type': 'svg',
            ''
            'action': 'download'
            }
    else: return no_update
"""
@app.callback(
    [Output('reads-slider', 'min'),
     Output('reads-slider', 'max'),
     Output('reads-slider', 'step'),
     Output('reads-slider', 'value')],
    Input('dropdown-update-dataset', 'value'))
def update_dataset(dataset_path):
    global initial_df, selected_node, fragments_sum
    initial_df = pd.read_csv(dataset_path)
    selected_node = None
    min_edge, max_edge = np.log(initial_df['nb_fragments'].min()), np.log(initial_df['nb_fragments'].max())
    fragments_sum = initial_df['nb_fragments'].sum()
    return min_edge - 0.01*min_edge, max_edge + 0.01*max_edge, max_edge/100, min_edge + 0.5*(max_edge-min_edge)

def open_browser():
    webbrowser.open_new_tab("http://localhost:8080")

if __name__ == '__main__':
    #open_browser();
    app.run_server(debug=True,port=8080,host='0.0.0.0');
