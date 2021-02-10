#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 10:19:40 2021

@author: malte
"""

from dash import Dash, no_update
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output

import os, json, webbrowser
import pandas as pd
import numpy as np

import dash_cytoscape as cyto
cyto.load_extra_layouts()

app = Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server

dir_path = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(dir_path, "assets", 'mystylesheet.json')) as json_file:
    stylesheet = json.load(json_file)
    
all_edges = pd.read_csv(os.path.join("~/Data/vibrio/coexpression/", 'edges_new.txt'))
all_nodes = pd.read_csv(os.path.join("~/Data/vibrio/coexpression/", 'nodes.txt'), sep="\t")
#all_edges = all_edges[:500000]

color_count = pd.DataFrame({'count':np.zeros(len(np.unique(all_nodes['color'])))}, index=np.unique(all_nodes['color']))
colors = {}
subgraphs = {}
for i, node in all_nodes.iterrows():
    colors[node['nodeName']] = node['color']
    color_count.loc[node['color']] += 1

color_count = color_count.sort_values(by=['count'])
color_order = color_count.index
for i, color in enumerate(color_order): subgraphs[color] = 'network {} ({})'.format(i+1, int(color_count.loc[color]))

max_corr = np.max(all_edges['weight'])

no_data = {"layout": {"xaxis": {"visible": False}, "yaxis": {"visible": False},
                      "annotations": [{"text": "No matching data found",
                                       "xref": "paper",
                                       "yref": "paper",
                                       "showarrow": False,
                                       "font": {"size": 28}}]}}

def cytoscape_data(gene, min_correlation, mode):
    
    jsondata = {"nodes": [], "edges": []}
    nodes = jsondata["nodes"]
    edges = jsondata["edges"]
    collected_nodes = set()
    
    if gene is None: return jsondata
    
    color = all_nodes.iloc[np.argmax(all_nodes['nodeName']==gene)]['color']
    corr_index = (all_edges['weight']>=min_correlation)
    if mode == 'subgraph':
        index = corr_index & ((all_edges['fromColor'] == color) & (all_edges['toColor'] == color))
        nodes_invloved = set(all_edges['toNode'][(all_edges['fromNode'] == gene) & corr_index]) | set(all_edges['fromNode'][(all_edges['toNode'] == gene) & corr_index])
    elif mode == 'all':
        index = corr_index & ((all_edges['fromNode'] == gene) | (all_edges['toNode'] == gene))
    elif mode == 'targets':
        index = corr_index & ((all_edges['fromColor'] == color) & (all_edges['toColor'] == color)) & ((all_edges['fromNode'] == gene) | (all_edges['toNode'] == gene))
        
    for i, coexpression in all_edges[index].iterrows():
        correlation = coexpression['weight']
        source = coexpression['fromNode']
        target = coexpression['toNode']
        if mode == 'subgraph':
            if (source in nodes_invloved and target == gene) or (target in nodes_invloved and source == gene): c = 'red'
            else: c='grey'
        else: c='grey'
        edges.append({'data': dict(source=source, target=target, correlation=correlation, color=c),
                      'classes': 'coexpression'}) 
        
        collected_nodes |= set([source, target])
            
    for node in collected_nodes:
        if mode == 'subgraph': 
            if node == gene: c = 'grey'
            else: c = colors[node]
        else: c = colors[node]
        nodes.append({'data':dict(id=node, name=node, color=c), 'classes':'gene'})
    
    return jsondata

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
                                html.H3(children="Coexpression Graph"),
                            ],
                        ),
                        html.Div(
                            className="container",
                            children=[
                                html.Div(
                                    className="control-block horizontal",
                                    children=[
                                        dcc.Dropdown(
                                            id='dropdown-update-layout',
                                            value='concentric',
                                            clearable=False,
                                            options=[
                                                {'label': name.capitalize(), 'value': name}
                                                for name in ['grid', 'random', 'circle', 'concentric']
                                            ]
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="control-block",
                                    children=[
                                        html.Div(
                                            className='control-element',
                                            children=[
                                                dcc.Dropdown(
                                                    placeholder="Search for a gene...",
                                                    value=None,
                                                    id='dropdown-gene',
                                                    clearable=False,
                                                    options=[
                                                        {'label':row['nodeName'], 'value':row['nodeName']} for i,row in all_nodes.iterrows()
                                                    ]
                                                ),
                                            ]
                                        ),
                                        html.Div(
                                            className='control-element',
                                            children=[
                                                #html.P("min reads:", className="line-label"),
                                                dcc.Slider(
                                                    id='slider-correlation',
                                                    min=0.1,
                                                    max=0.4,
                                                    step=(0.1)/100,
                                                    value=0.15,
                                                    marks={
                                                        0.1: '{}'.format(0.1),
                                                        0.2: '{}'.format(0.2),
                                                        0.3: '{}'.format(0.3),
                                                        0.4: '{}'.format(0.4)
                                                    }
                                                ),
                                            ]
                                        ),
                                        html.Div(
                                            className='control-element',
                                            children=[
                                                #html.P("min reads:", className="line-label"),
                                                dcc.RadioItems(
                                                    id='radio-filter',
                                                    options=[
                                                        {'label': 'all coexpressed', 'value': 'all'},
                                                        {'label': 'coexpressed in subnetwork', 'value': 'targets'},
                                                        {'label': 'complete subnetwork', 'value': 'subgraph'}
                                                    ],
                                                    value='all'
                                                )
                                            ]
                                        )
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
                            value='data', 
                            children=[
                                dcc.Tab(
                                    id="table-tab",
                                    className='custom-tab',
                                    selected_className='custom-tab--selected',
                                    label='Data',
                                    value='data',
                                    children=[
                                        html.Div(
                                            id="graph-container",
                                            className="container",
                                            children=[
                                                cyto.Cytoscape(
                                                    id='graph',
                                                    elements=cytoscape_data(None, 0.1, 'all'),
                                                    stylesheet=stylesheet,
                                                    responsive=True,
                                                    layout={'name':'concentric'}
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
                                                html.P("This is some explanation"),
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
    Output('graph', 'elements'),
    [Input('dropdown-gene', 'value'),
     Input('slider-correlation', 'value'),
     Input('radio-filter', 'value')])
def update_selected_data(gene, correlation, mode):
    return cytoscape_data(gene, correlation, mode)
    

@app.callback(
    Output('graph', 'layout'),
    Input('dropdown-update-layout', 'value'))
def update_layout(layout):
    return {
        'name': layout,
        'animate': False
    }

def open_browser():
    webbrowser.open_new_tab("http://localhost:8082")

if __name__ == '__main__':
    #open_browser();
    app.run_server(debug=False,port=8082,host='0.0.0.0');