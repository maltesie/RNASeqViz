from dash import Dash, html, dcc, dash_table, no_update
import dash_bio as dashbio
from dash.dependencies import Input, Output, State
import matplotlib.pyplot as plt
import matplotlib

import os
import pandas as pd
import networkx as nx
import json
import numpy as np
from itertools import combinations
import webbrowser

import dash_cytoscape as cyto
cyto.load_extra_layouts()

srna_types = ("sRNA", "ncRNA")
max_interactions = 300
fun_colors = ["DarkOliveGreen", "DarkRed", "DarkSlateBlue", "LightCoral", "Khaki"]

app = Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server

# Define Helper Functions

def filter_threshold(df, threshold):
    return df[np.log(df['nb_ints']) > threshold]

def filter_targets(df, node):
    return df[df['name1'].str.match(node) | df['name2'].str.match(node)]

def filter_cc(df, node):
    g = nx.Graph()
    g.add_edges_from(zip(df['name1'], df['name2']))
    for component in nx.connected_components(g):
        if node in component: break
    search_str = r"{}".format('|'.join(component))
    return df[df['name1'].str.match(search_str) | df['name2'].str.match(search_str)]

def filter_search(df, search_strings):
    index = np.zeros(len(df), dtype=bool)
    for search_string in search_strings:
        index |= df['name1'].str.contains(search_string) | df['name2'].str.contains(search_string)
    return df[index]

def table_data(df):
    return df[['type1', 'name1', 'name2', 'type2', 'in_libs', 'nb_ints']]

def fix_ig_label(label):
    if ':' in label:
        return '\n'.join(['IG']+label.split(':'))
    else: return label
    
def translate(value, leftMin, leftMax, rightMin, rightMax):
    leftSpan = leftMax - leftMin
    rightSpan = rightMax - rightMin
    valueScaled = min(max(float(value - leftMin) / float(leftSpan), -0.1), 1.1)
    return rightMin + (valueScaled * rightSpan)

def cb_edge_percentage(value):
    return int(round(translate(value, 0, 1, 17, 73)))

def cb_node_percentage(value):
    return int(round(translate(value, 0, 0.1, 14, 72)))

def type_by_function(n, t, fs):
    if fs is None or n not in multifun_data: return t
    s = "_".join([f for f in fs if f in multifun_data[n]])
    if len(s) == 0 : return t
    else: return s
    
def all_combinations(l):
    combs = []
    for i in range(1, len(l)+1):
        combs += combinations(l, i)
    return combs

def shorten(s):
    if len(s)>14:
        s=s[:14] + '...'
    return s

def cytoscape_data(df, norm, selected_functions):
    jsondata = {"nodes": [], "edges": []}
    nodes = jsondata["nodes"]
    edges = jsondata["edges"]
    fragment_count = {}
    current_norm = df['nb_ints'].sum()
    types = {}
    names = {}
    for i, interaction in df.iterrows():
        source, target, type_source, type_target, nb_ints, in_libs, pos1, pos2, min1, min2, max1, max2 = \
            interaction[['name1', 'name2', 'type1', 'type2', 'nb_ints', 'in_libs', 'relmean1', 'relmean2', 'relmin1', 'relmin2', 'relmax1', 'relmax2']]
        if type_source in srna_types: type_source = "sRNA"
        if type_target in srna_types: type_target = "sRNA"
        nb_ints /= norm
        pos, mi, ma = -1, -1, -1
        if ((type_target == "sRNA") != (type_source == "sRNA")): 
            edge_type = "srna_edge"
            if "IGR" in (type_target, type_source):
                edge_type = "other_edge"
            elif type_target in ("sRNA", "ncRNA"):
                pos, mi, ma = pos1, min1, max1
            else:
                pos, mi, ma = pos2, min2, max2     
            if pos == 0.0: pos=-0.1
            if pos == 1.0: pos=1.1
            if mi == 0.0: mi=-0.1
            if mi == 1.0: mi=1.1
            if ma == 0.0: ma=-0.1
            if ma == 1.0: ma = 1.1
        else: edge_type = "other_edge"
        
        source_id = source+'#'+type_source
        target_id = target+'#'+type_target
        edges.append({'data': dict(source=source_id, target=target_id, fragments=nb_ints, pos=pos, mi=mi, ma=ma, typ=edge_type, in_libs=in_libs),
                      'classes': edge_type}) 
        if source_id in fragment_count: 
            fragment_count[source_id] += nb_ints
        else: 
            fragment_count[source_id] = nb_ints
            types[source_id] = type_by_function(source, type_source, selected_functions)
            names[source_id] = source
        
        if target_id in fragment_count: 
            fragment_count[target_id] += nb_ints
        else: 
            fragment_count[target_id] = nb_ints
            types[target_id] = type_by_function(target, type_target, selected_functions)
            names[target_id] = target
        
        if i>max_interactions: break
            
    for node in fragment_count:
        nodes.append({'data':dict(id=node, name=fix_ig_label(names[node]), fragments=fragment_count[node], current_fragments=fragment_count[node]*norm/current_norm, typ=types[node]), 'classes':types[node]})
    
    return jsondata

def circos_data(df):
    data = []
    m = df['nb_ints'].max()
    cmap = plt.get_cmap('binary')
    for i, interaction in df.iterrows():
        source, target, nb_ints, start1, stop1, start2, stop2, chr1, chr2 = \
            interaction[['name1', 'name2', 'nb_ints', 'meanleft1', 'meanright1', 'meanleft2', 'meanright2', 'ref1', 'ref2']]
        rgb = matplotlib.colors.rgb2hex(cmap((nb_ints/m + 3)/4))
        data.append({
            'color': rgb,
            'nb_ints': nb_ints,
            "source": {
                "id": chr1,
                "start": abs(start1),
                "end": abs(start1)+1500
            },
            "target": {
                "id": chr2,
                "start": abs(start2),
                "end": abs(start2)+1500
            }
        })
    return data

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
min_edge, max_edge = np.log(initial_df['nb_ints'].min()), np.log(initial_df['nb_ints'].max())
fragments_sum = initial_df['nb_ints'].sum()
elements = cytoscape_data(initial_df[:100], fragments_sum, None)

genes = [{"label":v, "value":v} for v in np.sort(np.unique(np.vstack((initial_df["name1"],initial_df["name2"]))))]

multifun_trans = {row["multi_fun"].replace(".", "-"):row["common_name"] for i, row in pd.read_csv(os.path.join(dir_path, "assets", "multifun_categories.txt"), sep='\t').iterrows()}
multifun_data = {row["name"]:"_".join([mf.replace(".", "-") for mf in row["multi_fun"].split(":")]) for i, row in pd.read_csv(os.path.join(dir_path, "assets", "multifun_genes.txt"), sep='\t').iterrows() if row["multi_fun"].startswith("BC-")}
multifun_items = [{"value":key, "label":"{}: {}".format(key, val)} for key, val in multifun_trans.items()]

selected_fun_genes = []

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
                                                        for name in ['cose', 'grid', 'random', 'circle', 'concentric']
                                                    ]
                                                )
                                            ]
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="controls-block",
                                    children=[
                                        html.P("Color by function:"),
                                        dcc.Dropdown(
                                            id='function-multi-select',
                                            options=[{"label":"top 5 in current graph", "value":"all"}]+multifun_items,
                                            multi=True
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
                                                {'label': 'restrict to selected functions', 'value': 'fun'},
                                                {'label': 'only targets of selected', 'value': 'targets'},
                                                {'label': 'connected component of selected', 'value': 'cc'}
                                            ],
                                            value='all'
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="controls-block",
                                    children=[
                                        html.P("Search:"),
                                        dcc.Dropdown(
                                            id="gene-multi-select",
                                            options=genes,
                                            multi=True
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
                                                    layout={'name':'random'},
                                                    minZoom=0.1,
                                                    maxZoom=5.0
                                                ),
                                                html.Div(
                                                    children=[html.Button('DOWNLOAD SVG', id='save-svg', n_clicks=0), html.Div(className="spacer"), html.Div(id="legend-container"), html.Div(className="spacer")])
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
                                                        'gap' : 0.01,
                                                        'cornerRadius':5,
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
                                                          "id": "NC_000913.3",
                                                          "label": "",
                                                          "color": "#999999",
                                                          "len": 4641652
                                                        }],
                                                    selectEvent={"0": "hover", "1": "click", "2": "both"},
                                                    tracks=[{
                                                        'type': 'CHORDS',
                                                        'data': circos_data(initial_df),
                                                        'config': {
                                                            'opacity': 0.9,
                                                            'color': {'name': 'color'},
                                                            'tooltipContent': {'name':'nb_ints'}
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
                                                ),
                                                html.Div(
                                                    [
                                                        html.Button("DOWNLOAD CSV", id="btn_csv"),
                                                        dcc.Download(id="download-dataframe-csv"),
                                                    ]
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
                                                    id = "test-block",
                                                    children=[
                                                        html.P("", className="text-block")
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
    Output("download-dataframe-csv", "data"),
    Input("btn_csv", "n_clicks"),
    prevent_initial_call=True,
)
def func(n_clicks):
    return dcc.send_data_frame(filtered_df.to_csv, "interactions.csv")

@app.callback(
    [Output("graph", "stylesheet"),
    Output('table', 'data'),
    Output('graph', 'elements'),
    Output('my-dashbio-circos', 'tracks'),
    Output('reads-slider', 'marks'),
    Output('legend-container', 'children'),
    Output('function-multi-select', 'value'),
    Output('test-block', 'children')],
    [Input('reads-slider', 'value'),
    Input('filter-radio', 'value'),
    Input('gene-multi-select', 'value'),
    Input('function-multi-select', 'value')])
def update_selected_data(slider_value, radio_filter, search_strings, functions):
    global filtered_df, selected_fun_genes
    fun_output = no_update
    if functions is None: selected_fun_genes = []
    else:
        if "all" in functions: 
            all_names = np.unique(np.hstack((filtered_df["name1"][filtered_df["type1"]=="CDS"], filtered_df["name2"][filtered_df["type2"]=="CDS"])))
            fun_output = []
            count = {c:0 for c in multifun_trans}
            for name in all_names: 
                if name in multifun_data: 
                    for mf in multifun_data[name].split("_"): count[mf] += 1
            functions = fun_output = sorted(count, key=count.get, reverse=True)[:5]    
        selected_fun_genes = [g for g,mfs in multifun_data.items() if any([f in mfs for f in functions])]
    filtered_df = filter_threshold(initial_df, slider_value)
    if (radio_filter == 'targets') and (selected_node is not None): filtered_df = filter_targets(filtered_df, selected_node)
    elif (radio_filter == 'cc') and (selected_node is not None): filtered_df = filter_cc(filtered_df, selected_node)
    elif radio_filter == 'fun': filtered_df = filter_search(filtered_df, selected_fun_genes)
    if search_strings is not None and len(search_strings) > 0: filtered_df = filter_search(filtered_df, search_strings)
    tracks=[{
        'type': 'CHORDS',
        'data': circos_data(filtered_df),
        'config': {
            'logScale': False,
            'opacity': 0.9,
            'color': {'name': 'color'},
            'tooltipContent': {'name': 'nb_ints'}
        }
    }]
    if functions is not None:
        my_stylesheet = stylesheet + [
            {"selector":"."+"_".join(fun_combination), 
             "style":{
                 "background-fill": "linear-gradient",
                 "background-gradient-stop-colors": " ".join(color_combination),
                 "background-gradient-direction": "to-right",
                 "background-blacken":"-0.2",
                 "font-weight": "bold"
                 }
             } 
            for fun_combination, color_combination in zip(all_combinations(functions), all_combinations(fun_colors[0:len(functions)]))]
    else:
        my_stylesheet = stylesheet
    
    if functions is not None: legend = [html.Div(className="horizontal", children=[html.Div(className="color-legend-item", 
                                                                                            style={"background-color":c}), 
                                                                                   html.P(f, className="color-legend-text")]) for c, f in zip(fun_colors[:len(functions)], functions)]
    else: legend = []
    
    return my_stylesheet, table_data(filtered_df).to_dict('records'), \
            cytoscape_data(filtered_df, fragments_sum, functions), \
            tracks, {slider_value: '{} ({})'.format(int(np.exp(slider_value)), len(filtered_df))}, \
            legend, fun_output, [html.P(str(fun_output)+str(functions))]

@app.callback(
    [Output('filter-radio', 'options'),
     Output('test-output', 'children'),
     Output('filter-radio', 'value')],
    [Input('graph', 'selectedNodeData'),
     Input('graph', 'selectedEdgeData')],
    State('filter-radio', 'value'))
def set_selected_element(node_data, edge_data, radio_value):
    global selected_node
    if (node_data is None) or (node_data == []):
        selected_node = None
        radio_return = [{'label': 'all interactions', 'value': 'all'},
                        {'label': 'restrict to selected functions', 'value': 'fun'},
                        {'label': 'only interactions with selected RNA', 'value': 'targets'},
                        {'label': 'connected component of selected RNA', 'value': 'cc'}]
        text_return = ["Select a node or edge in the graph."]
        if radio_value != 'all': value_return = 'all'
        else: value_return = no_update
        if (edge_data is None) or (edge_data == []):
            text_return = ["Select an interaction to display further information."]
        else:  
            edge_data = edge_data[0]
            text_return = ["{} -> {}".format(edge_data['source'].split("#")[0], edge_data['target'].split("#")[0]), html.Br(), 
                           "# of reads: {}".format(int(fragments_sum*float(edge_data['fragments'])))]
            value_return = no_update
            if edge_data['typ'] == 'srna_edge':
                text_return.append(html.Div(className="horizontal deflate", children=[html.P("5'UTR", className="cb-left"),
                                                                                        html.Div(
                                                                                            id = "colorbar-edges",
                                                                                            className="controls-block",
                                                                                            children=[
                                                                                                html.P("CDS", className="cb-middle"),
                                                                                                html.Div(
                                                                                                    className='colorbar-arrow',
                                                                                                    style={
                                                                                                        'left': '{}%'.format(cb_edge_percentage(edge_data['mi'])),
                                                                                                    }
                                                                                                ),
                                                                                                html.Div(
                                                                                                    className='colorbar-arrow',
                                                                                                    style={
                                                                                                        'left': '{}%'.format(cb_edge_percentage(edge_data['ma'])),
                                                                                                    }
                                                                                                ),
                                                                                                html.Div(
                                                                                                    className='colorbar-arrow mean',
                                                                                                    style={
                                                                                                        'left': '{}%'.format(cb_edge_percentage(edge_data['pos']))
                                                                                                    }
                                                                                                )
                                                                                            ]
                                                                                        ),
                                                                                        html.P("3'UTR", className="cb-right")]))

    else: 
        node_data = node_data[0]
        selected_node = node_data["name"]
        radio_return = [{'label': 'all interactions', 'value': 'all'},
                        {'label': 'restrict to selected functions', 'value': 'fun'},
                        {'label': 'only interactions with {}'.format(selected_node), 'value': 'targets'},
                        {'label': 'connected component of {}'.format(selected_node), 'value': 'cc'}]
        selected_node_interactions = filter_targets(filtered_df, selected_node)
        nb_targets = len(set(list(selected_node_interactions['name1'])+list(selected_node_interactions['name2']))) - 1
        text_return = ["{} is involved in {} interactions with {} partners.".format(selected_node, len(selected_node_interactions), nb_targets)]
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
                                                                                                    'left': '{}%'.format(cb_node_percentage(selected_node_interactions['nb_ints'].sum()/filtered_df['nb_ints'].sum())),
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
            'componentSpacing': 10,
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

@app.callback(
    [Output('reads-slider', 'min'),
     Output('reads-slider', 'max'),
     Output('reads-slider', 'step'),
     Output('reads-slider', 'value')],
    Input('dropdown-update-dataset', 'value'))
def update_dataset(dataset_path):
    global initial_df, selected_node, fragments_sum, genes
    initial_df = pd.read_csv(dataset_path)
    genes = [{"label":v, "value":v} for v in np.sort(np.unique(np.vstack((initial_df["name1"],initial_df["name2"]))))]
    selected_node = None
    min_edge, max_edge = np.log(initial_df['nb_ints'].min()), np.log(initial_df['nb_ints'].max())
    fragments_sum = initial_df['nb_ints'].sum()
    return min_edge - 0.01*min_edge, max_edge + 0.01*max_edge, max_edge/100, min_edge + 0.5*(max_edge-min_edge)

def open_browser():
    webbrowser.open_new_tab("http://localhost:8080")

if __name__ == '__main__':
    #open_browser();
    app.run_server(debug=True,port=8081,host='0.0.0.0');
