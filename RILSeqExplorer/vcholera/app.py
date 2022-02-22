from dash import Dash, html, dcc, dash_table, no_update, callback_context
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

functional_annotation_db = "KEGG Orthology" # "Multifun"
srna_types = ("sRNA", "ncRNA")
fun_colors = ["DarkOliveGreen", "DarkRed", "DarkSlateBlue", "LightCoral", "Khaki", "Blue", "Chartreuse", "Cyan", "BlueViolet", "DarkMagenta"]

app = Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server

# Define Helper Functions

def hierarchical_sort(sorted_items):
    hier_sorted = []
    while sorted_items:
        s1 = sorted_items.pop(0)
        hier_sorted.append(s1)
        while True:
            s2 = s1
            l = len(sorted_items)
            for i, other_s in enumerate(sorted_items):
                if other_s.startswith(s2):
                    if s2 == s1:
                        s2 = other_s
                    hier_sorted.append(other_s)
                    sorted_items.remove(other_s)
            if len(sorted_items) == l: break
    return hier_sorted

len_trans = {5:"", 8:"- ", 11:"- - "}

def filter_df(df, search_strings, max_interactions, slider_value, functions, checklist):
    filtered_df = df.copy()
    if search_strings is not None and len(search_strings) > 0: filtered_df = filter_search(filtered_df, search_strings)
    filtered_df = filter_threshold(filtered_df, slider_value)[:max_interactions]

    if checklist is None: checklist = []
    
    all_names = np.unique(np.hstack((filtered_df["name1"][filtered_df["type1"].isin(["CDS", "5UTR", "3UTR"])], filtered_df["name2"][filtered_df["type2"].isin(["CDS", "5UTR", "3UTR"])])))
    #if top10 or top10g:
    count = {c:0 for c in multifun_trans}
    for name in all_names: 
        if name in multifun_data: 
            for mf in multifun_data[name].split("_"): 
                count[mf] += 1
                #_, a, b, c = mf.split("-")
                #count["KO" + "-" + a + "-" + b] += 1
                #count["KO" + "-" + a] += 1
        
    sorted_items = sorted({k:v for k,v in count.items() if v>0}, key=count.get, reverse=True)
    
    select_options = [{"label":"[{}] {}".format(count[a], multifun_trans[a]), "value":a}  for a in sorted_items]
    #if top10g: functions = [f for i,f in enumerate(functions) if not any(((count[ff] > count[f]) and (ff+'-' in f)) or ((count[ff] == count[f]) and (f+'-' in ff)) for ff in functions)]
    #functions = functions[:10]
    
    selected_fun_genes = []
    fun2color = dict()
    if functions is not None: 
        fun2color = {f:c for f,c in zip(functions, fun_colors)}
        selected_fun_genes = [n for n in all_names if ((n in multifun_data) and any(any(ff == f or f+'-' in ff for ff in multifun_data[n].split("_")) for f in functions))] 
        if ("restrict_fun" in checklist): filtered_df = filter_search(filtered_df, selected_fun_genes)
        
    return filtered_df, functions, fun2color, select_options, selected_fun_genes

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
        index |= df['name1'].str.match(search_string) | df['name2'].str.match(search_string)
    return df[index]

def table_data(df):
    return df[['type1', 'name1', 'name2', 'type2', 'in_libs', 'nb_ints']]

def fix_ig_label(label):
    if ':' in label:
        return '\n'.join(['IG']+label.split(':'))
    #elif label.startswith('VC'): 
    #    if label[2] == 'A': return 'VCA\n'+label[3:]
    #    else: return 'VC\n'+label[2:]
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
    
def all_combinations(l):
    combs = []
    for i in range(1, len(l)+1):
        combs += combinations(l, i)
    return combs

def shorten(s):
    if len(s)>14:
        s=s[:14] + '...'
    return s

def break_lines(s, nchars=18):
    a = s.split(" ")
    return_s = ""
    c = 0
    for w in a:
        c += len(w)
        if c > nchars:
            c = len(w)
            return_s += "\n"
        return_s += w
        if c <= nchars:
            c += 1
            return_s += " "
    return return_s

def type_by_function2(n, t, fs):
    if fs is None or n not in multifun_data or not any(f in multifun_data[n] for f in fs): return t
    return multifun_data[n]
    
def type_by_function(n, t):
    if n not in multifun_data: return t
    return multifun_data[n]

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
            elif type_target in srna_types:
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
        
        source_id = source
        target_id = target
        edges.append({'data': dict(id=target+source, source=source_id, target=target_id, fragments=nb_ints, norm_fragments=nb_ints*norm/current_norm, pos=pos, mi=mi, ma=ma, typ=edge_type, in_libs=in_libs),
                      'classes': edge_type}) 
        if source_id in fragment_count: 
            fragment_count[source_id] += nb_ints
        else: 
            fragment_count[source_id] = nb_ints
            types[source_id] = type_by_function2(source, type_source, selected_functions)
            names[source_id] = source
        
        if target_id in fragment_count: 
            fragment_count[target_id] += nb_ints
        else: 
            fragment_count[target_id] = nb_ints
            types[target_id] = type_by_function2(target, type_target, selected_functions)
            names[target_id] = target
            
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
initial_dfs = []
dir_path = os.path.dirname(os.path.realpath(__file__))
csv_index = 0
csv_trans = {}
fragments_sums = []
for file in os.listdir(os.path.join(dir_path, "assets")):
    if file.endswith(".csv"):
        csv_path = os.path.abspath(os.path.join(dir_path, "assets", file))
        dataset_paths.append(csv_path)
        initial_dfs.append(pd.read_csv(csv_path))
        fragments_sums.append(initial_dfs[-1]['nb_ints'].sum())
        csv_trans[csv_path] = csv_index
        csv_index += 1
        
        
with open(os.path.join(dir_path, "assets", 'mystylesheet.json')) as json_file:
    stylesheet = json.load(json_file)

if functional_annotation_db == "KEGG Orthology":
    gene_name = "name"
    identifier = "hierarchy"
    description = "description"
    gene_categories = "ko_category"
    fname_genes = "ko_genes.txt"
    fname_categories = "ko_categories.txt"
    category_prefix = "KO-"
elif functional_annotation_db == "Multifun":
    gene_name = "name"
    identifier = "multi_fun"
    description = "common_name"
    gene_categories = "multi_fun"
    fname_genes = "multifun_genes.txt"
    fname_categories = "multifun_categories.txt"
    category_prefix = "BC-"

multifun_trans = {row[identifier].replace(".", "-"):row[description] for i, row in pd.read_csv(os.path.join(dir_path, "assets", fname_categories), sep='\t').iterrows()}
multifun_data = {row[gene_name]:"_".join([mf.replace(".", "-") for mf in row[gene_categories].split(":")]) for i, row in pd.read_csv(os.path.join(dir_path, "assets", fname_genes), sep='\t').iterrows() if row[gene_categories].startswith(category_prefix)}
multifun_set = np.unique(np.hstack([[mf.replace(".", "-") for mf in row[gene_categories].split(":")] for i, row in pd.read_csv(os.path.join(dir_path, "assets", fname_genes), sep='\t').iterrows() if row[gene_categories].startswith(category_prefix)]))
multifun_items = [{"value":key, "label":"[0] "+val} for key, val in multifun_trans.items() if key in multifun_set]

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
                                        html.P("Search:"),
                                        dcc.Dropdown(
                                            id="gene-multi-select",
                                            options=[],
                                            multi=True
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="controls-block",
                                    children=[
                                        html.Div(
                                            className="horizontal",
                                            children=[
                                                html.P("Minimum # of reads:"),
                                                dcc.Slider(
                                                    id='reads-slider',
                                                    min=0,
                                                    max=100,
                                                    step=1,
                                                    value=50
                                                ),
                                            ]
                                        ),
                                        html.Div(
                                            className="horizontal",
                                            style={"padding-top":"10px"},
                                            children=[
                                                html.P(id="nb-interactions-text", style={"padding-right":"5px", "padding-top":"2px"}),
                                                dcc.Input(id="max-interactions", type="number", value=300, style={"max-height":"18px", "min-width":"80px"})
                                            ]
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="controls-block",
                                    children=[
                                        html.P("Annotated Functions:"),
                                        dcc.Dropdown(
                                            id='function-multi-select',
                                            options=multifun_items,
                                            clearable=False,
                                            multi=True
                                        ),
                                        dcc.Checklist(
                                            id="color-checklist",
                                            options=[
                                                {"label":"Restrict current graph to selected functions", "value":"restrict_fun"}
                                            ],
                                            style={"padding-top":"12px", "display":"flex", "flex-direction":"column"}
                                        )
                                    ]
                                ),
                                
                                html.Div(
                                    className="controls-block",
                                    children=[
                                        html.P(children=["No interaction selected."], id="info-output")
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
                                                    elements=[],
                                                    stylesheet=stylesheet,
                                                    responsive=True,
                                                    layout={'name':'random'},
                                                    minZoom=0.1,
                                                    maxZoom=5.0
                                                ),
                                                html.Div(
                                                    children=[html.Button('DOWNLOAD SVG', id='save-svg', n_clicks=0), html.Div(style={"height":"8%"}), html.Div(id="legend-container"), html.Div(style={"height":"100%"})])
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
                                                        'gap' : 0.003,
                                                        'cornerRadius':5,
                                                        'innerRadius':290,
                                                        'ticks': {
                                                            'display': True, 
                                                            'spacing': 100000,
                                                            'color': '#000',
                                                            'labelDenominator': 1000000,
                                                            'labelSuffix': ' Mb',
                                                            'labelFont': 'Arial',
                                                            'majorSpacing':5,
                                                            'minorSpacing':1,
                                                            'labelSpacing':5
                                                            }
                                                        },
                                                    layout=[
                                                        {
                                                          "id": "NC_002505",
                                                          "label": "",
                                                          "color": "#009933",
                                                          "len": 2961149
                                                        },
                                                        {
                                                          "id": "NC_002506",
                                                          "label": "",
                                                          "color": "#99cc33",
                                                          "len": 1072315
                                                        }
                                                    ],
                                                    selectEvent={"0": "hover", "1": "click", "2": "both"},
                                                    tracks=[{
                                                        'type': 'CHORDS',
                                                        'data': [],
                                                        'config': {
                                                            'opacity': 0.9,
                                                            'color': {'name': 'color'},
                                                            'tooltipContent': {'name':'nb_ints'}
                                                        }
                                                    }]
                                                )
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
                                                    columns=[{"name": i, "id": i} for i in ["name1", "type1", "name2", "type2", "nb_ints", "in_libs"]],
                                                    style_cell={
                                                        'height': 'auto',
                                                        'minWidth': '140px', 'width': '140px', 'maxWidth': '140px',
                                                        'whiteSpace': 'normal'
                                                    }
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
    [State('gene-multi-select', 'value'),
     State('max-interactions', 'value'),
     State('reads-slider', 'value'),
     State('function-multi-select', 'value'),
     State('color-checklist', 'value'),
     State('dropdown-update-dataset', 'value')],
    prevent_initial_call=True
)
def func(n_clicks, search_strings, max_interactions, slider_value, functions, checklist, dataset):
    return dcc.send_data_frame(filter_df(initial_dfs[csv_trans[dataset]], search_strings, max_interactions, slider_value, functions, checklist)[0].to_csv, "interactions.csv")

@app.callback(
    [Output('graph', 'layout'),
    Output('graph', 'elements'),
    Output("graph", 'stylesheet'),
    Output('table', 'data'),
    Output('my-dashbio-circos', 'tracks'),
    Output('reads-slider', 'marks'),
    Output('nb-interactions-text', 'children'),
    Output('function-multi-select', 'options'),
    Output('legend-container', 'children')],
    [Input('dropdown-update-layout', 'value'),
    Input('reads-slider', 'value'),
    Input('gene-multi-select', 'value'),
    Input('color-checklist', 'value'),
    Input('max-interactions', 'value'),
    Input('dropdown-update-dataset', 'value'),
    Input('function-multi-select', 'value')],
    State('graph', 'elements'),
    prevent_initial_call=True
    )
def update_selected_data(layout_value, slider_value, search_strings, checklist, max_interactions, dataset, functions):
    
    filtered_df, functions, fun2color, fun_items, selected_fun_genes = filter_df(initial_dfs[csv_trans[dataset]], search_strings, max_interactions, slider_value, functions, checklist)

    all_unique_fun_strings = np.unique([multifun_data[n] for n in selected_fun_genes])
    my_stylesheet = [
        {"selector":"."+fun_combination, 
         "style":{
             "background-fill": "linear-gradient",
             "background-gradient-stop-colors": " ".join([fun2color[fun] for fun in fun2color if any(fun==f for f in fun_combination.split("_"))]),
             "background-gradient-direction": "to-right",
             "background-blacken":"-0.2",
             "font-weight": "bold"
             }
         } 
        for fun_combination in all_unique_fun_strings]

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
    
    if functions is not None: legend = [html.Div(className="horizontal centered", children=[html.Div(className="color-legend-item", style={"background-color":c}), 
                                                                                   html.P(break_lines(multifun_trans[f]), className="color-legend-text")]) for c, f in zip(fun_colors, functions)]
    else: legend = []
    
    table = table_data(filtered_df).to_dict('records')
    graph = cytoscape_data(filtered_df, fragments_sums[csv_trans[dataset]], functions)
    circos = tracks
    slider_value = {slider_value: '{}'.format(int(round(np.exp(slider_value)))+1)}
    slider_text = ["Current # of interactions: {} / ".format(len(filtered_df))]

    if layout_value == 'cose-bilkent':
        layout = {
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
    elif layout_value == 'concentric':
        layout = {
            'name':'concentric',
            'animate':False
        }
    else:
        layout = {
            'name': layout_value,
            'animate': False
        }

    if callback_context.triggered:
        for t in callback_context.triggered:
            if t["prop_id"] == "function-multi-select.value":
                layout = {
                    'name': 'preset',
                    'animate': False
                }### LOOK INTO THIS
                break
    #if graph_elements: legend = [html.P(str(graph_elements))]
    return layout, graph, stylesheet+my_stylesheet, table, circos, slider_value, slider_text, fun_items, legend
    
@app.callback(
    Output('info-output', 'children'),
    [Input('graph', 'selectedNodeData'),
     Input('graph', 'selectedEdgeData')],
    [State('gene-multi-select', 'value'),
     State('max-interactions', 'value'),
     State('reads-slider', 'value'),
     State('function-multi-select', 'value'),
     State('color-checklist', 'value'),
     State('dropdown-update-dataset', 'value')]
    )
def set_selected_element(node_data, edge_data, search_strings, max_interactions, slider_value, functions, checklist, dataset):
    global selected_node
    if (node_data is None) or (node_data == []):
        selected_node = None
        text_return = ["Select a node or edge in the graph."]
        if (edge_data is None) or (edge_data == []):
            text_return = ["Select an interaction (edge) or target (node) to display details. The size of the nodes corresponds to their share of interactions in the current graph. The color of the edges indicates the average binding position from red for 5'UTR to blue for 3'UTR."]
        else:  
            edge_data = edge_data[0]
            text_return = ["{} -> {}".format(edge_data['source'].split("#")[0], edge_data['target'].split("#")[0]), html.Br(), 
                           "# of reads: {}".format(int(fragments_sums[csv_trans[dataset]]*float(edge_data['fragments'])))]
            if edge_data['typ'] == 'srna_edge':
                text_return += [html.Br(), html.Br(), "Average chimera position in target:", html.Br()]
                text_return.append(html.Div(className="horizontal deflate", children=[html.P("5'UTR", className="cb-left"),
                                                                                        html.Div(
                                                                                            id = "colorbar-edges",
                                                                                            className="controls-block",
                                                                                            children=[
                                                                                                html.P("CDS", className="cb-middle"),
                                                                                                
                                                                                                html.Div(
                                                                                                    className='colorbar-arrow',
                                                                                                    style={
                                                                                                        'left': '{}%'.format(cb_edge_percentage(edge_data['pos']))
                                                                                                    }
                                                                                                )
                                                                                            ]
                                                                                        ),
                                                                                        html.P("3'UTR", className="cb-right")]))

    else: 
        node_data = node_data[0]
        selected_node = node_data["id"]
        selected_node_interactions = filter_targets(filter_df(initial_dfs[csv_trans[dataset]], search_strings, max_interactions, slider_value, functions, checklist)[0], selected_node)
        nb_targets = len(set(list(selected_node_interactions['name1'])+list(selected_node_interactions['name2']))) - 1
        if selected_node in multifun_data: node_funs = multifun_data[selected_node].split("_")
        else: node_funs= []
        t = ""
        if nb_targets > 1: t = "s"
        text_return = ["{} is involved in {} interaction{}.".format(selected_node, nb_targets, t)]
        if len(node_funs) > 0: text_return += [html.Br(), html.Br(), "{} classes:".format(functional_annotation_db)]
        else: text_return += [html.Br(), html.Br(), "No {} classes are annotated.".format(functional_annotation_db)]
        for nf in node_funs:
            text_return += [html.Br(), multifun_trans[nf]]
        
    return text_return
    
@app.callback(
    Output("graph", "generateImage"),
    Input('save-svg', 'n_clicks'))
def get_image(clicks):
    if clicks>0:
        return {
            'type': 'svg',
            'action': 'download'
            }
    else: return no_update

'''
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
'''
@app.callback(
    [Output('reads-slider', 'min'),
     Output('reads-slider', 'max'),
     Output('reads-slider', 'step'),
     Output('reads-slider', 'value'),
     Output('gene-multi-select', 'options')],
    Input('dropdown-update-dataset', 'value'))
def update_dataset(dataset_path):
    genes = [{'label':v, 'value':v} for v in np.sort(np.unique(np.vstack((initial_dfs[csv_trans[dataset_path]]["name1"],initial_dfs[csv_trans[dataset_path]]["name2"]))))]
    min_edge, max_edge = np.log(initial_dfs[csv_trans[dataset_path]]['nb_ints'].min()), np.log(initial_dfs[csv_trans[dataset_path]]['nb_ints'].max())
    return min_edge - 0.01*min_edge, max_edge + 0.01*max_edge, max_edge/100, min_edge + 0.5*(max_edge-min_edge), genes

def open_browser():
    webbrowser.open_new_tab("http://localhost:8080")

if __name__ == '__main__':
    #open_browser();
    app.run_server(debug=False,port=8081,host='0.0.0.0');
