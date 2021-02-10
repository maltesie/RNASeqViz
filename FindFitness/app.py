#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 10:19:40 2021

@author: malte
"""

from dash import Dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output

import os
import pandas as pd
import numpy as np
import plotly.express as px
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

no_data = {"layout": {"xaxis": {"visible": False}, "yaxis": {"visible": False},
                      "annotations": [{"text": "No matching data found",
                                       "xref": "paper",
                                       "yref": "paper",
                                       "showarrow": False,
                                       "font": {"size": 28}}]}}

all_data = {'caulo':{}, 'ecoli':{}}
for species in all_data:
    data = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assets', species, 'strain_fit.tab'), sep="\t")
    meta = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assets', species, 'fit_quality.tab'), sep="\t", index_col='name')
    meta = meta[meta['u']]
    meta['group'] = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assets', species, 'expsUsed'), sep="\t", index_col='name')['Group']
    
    data = data[['barcode','rcbarcode','scaffold','strand','pos','locusId','f','enoughT0'] + list(meta.index)]
    #data = data[data['f'].isna()]
    
    #if species == 'caulo':
    #    print(list(meta.index))
    #    print(data[['set1IT013', 'set1IT030', 'set1IT047']][0:100])
    
    #if species == 'ecoli':
    #    for i,row in data[['pos', 'f', 'locusId']][((4230000 <= np.abs(data['pos'])) & (np.abs(data['pos']) <= 4235000))].iterrows():
    #        print(row['pos'], " ", row['f'], " ", row['locusId'])
    
    groups = {}
    for group in np.unique(meta['group']):
        experiments = meta[meta['group'] == group].index
        data[group] = data[experiments].mean(axis=1)
        groups[group] = list(experiments)
    
    all_data[species]['data'] = data.copy()
    all_data[species]['meta'] = meta.copy()
    all_data[species]['groups'] = groups.copy()
    
annotation = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assets', 'caulo', 'NC_011916.gff'), sep="\t", skiprows=5)
annotation['locus_tag'] = np.roll([row['aux'].split(';')[-1][10:] for i, row in annotation.iterrows()], 1)
annotation = annotation[annotation['type']=='ncRNA']
annotation.reset_index()
#annotation.loc[annotation['strand']=='-','start'], annotation.loc[annotation['strand']=='-','stop'] = -1 * annotation[annotation['strand']=='-']['stop'], -1 * annotation[annotation['strand']=='-']['start']
all_data['caulo']['annotation'] = annotation.copy()

annotation = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assets', 'ecoli', 'srna_mg1655.csv'), sep=",")
#annotation.loc[annotation['strand']=='-','start'], annotation.loc[annotation['strand']=='-','stop'] = -1 * annotation[annotation['strand']=='-']['stop'], -1 * annotation[annotation['strand']=='-']['start']
#print(annotation['start'], annotation['stop'])
all_data['ecoli']['annotation'] = annotation.copy()

description = {**{'carbon source': 'all carbon source', 'nitrogen source': 'all nitrogen source', 'stress': 'all stress', 'pH':'all pH', 'pye':'all pye'},
               **{i:row['short'] for i, row in all_data['caulo']['meta'].iterrows()}}
all_data['caulo']['description'] = description.copy()

description = {**{'carbon source': 'all carbon source', 'nitrogen source': 'all nitrogen source', 'stress': 'all stress', 'pH':'all pH', 'pye':'all pye', 'sulfur source': 'all sulfur source', 'motility':'all motility'},
               **{i:row['short'] for i, row in all_data['ecoli']['meta'].iterrows()}}
all_data['ecoli']['description'] = description.copy()

#print(all_data['ecoli']['data'][(all_data['ecoli']['data']['pos'] <= -1995393) & (all_data['ecoli']['data']['pos'] >= -1995993)])

def start_end_from(name):
    index = np.argmax(annotation['locus_tag']==name)
    return annotation.iloc[index]['start'], annotation.iloc[index]['stop']


def data_plot(experiment_names, srna_name):
    locus_start, locus_end = start_end_from(srna_name)
    l = locus_end - locus_start
    index = (locus_start <= data['pos']) & (data['pos'] <= locus_end)
    rel_pos = []
    fitness = []
    condition = []
    belongs_to_exp = []
    for name in experiment_names:
        if name in groups:
            rel_pos += list((data[index]['pos'] - locus_start) / l) * len(groups[name])
            fitness += list(data[index][groups[name]].values.flatten())
            condition += [description[name]] * index.sum() * len(groups[name])
            belongs_to_exp += [description[n] for n in groups[name]] * index.sum()
        else:
            rel_pos += list((data[index]['pos'] - locus_start) / l)
            fitness += list(data[index][name])
            condition += [description[name]] * index.sum()
            belongs_to_exp += [description[name]] * index.sum()
    df = pd.DataFrame({'relative position':rel_pos, 'fitness':fitness, 'condition':condition, 'exp':belongs_to_exp})
    if len(df.index) > 0 : 
        fig = px.scatter(df, x='relative position', y='fitness', color='condition', hover_data=['exp'], title=srna_name, trendline="ols")
        fig.update_layout(clickmode='event+select')
        fig.update_xaxes(tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1.0], range=[0,1])
        fig.update_traces(marker_size=10)
        return fig
    else:
        return no_data

def top_ten_plot(experiment_names, mode='worst'):
    exps = set()
    title = []
    for exp in experiment_names:
        if exp in groups:
            exps |= set(groups[exp])
        else:
            exps |= set([exp])
        title.append(description[exp])
    title = ', '.join(title)
    exps=list(exps)
    avgs = []
    for srna_name in annotation['locus_tag']:
        locus_start, locus_end = start_end_from(srna_name)
        index = (locus_start <= data['pos']) & (data['pos'] <= locus_end)
        if index.sum() > 0: 
            avgs.append(data[index][exps].values.mean())
        else: avgs.append(0)
    if mode == 'top': index = np.argsort(avgs)[-5:][::-1]
    elif mode == 'worst': index = np.argsort(avgs)[:5]
    else: raise AssertionError('wrong mode')
    srna_names = annotation.iloc[index]['locus_tag'].values
    rel_pos = []
    fitness = []
    names = []
    belongs_to_exp = []
    for srna_name in srna_names:
        locus_start, locus_end = start_end_from(srna_name)
        l = locus_end - locus_start
        index = (locus_start <= data['pos']) & (data['pos'] <= locus_end)
        rel_pos += list(((data[index]['pos'] - locus_start) / l).values.flatten().repeat(len(exps)))
        fitness += list(data[index][exps].values.flatten())
        names += [srna_name] * index.sum() * len(exps)
        belongs_to_exp += [description[n] for n in exps] * index.sum()
    df = pd.DataFrame({'relative position':rel_pos, 'fitness':fitness, 'srna':names, 'condition':belongs_to_exp})
    if len(df.index) > 0: 
        fig = px.scatter(df, x='relative position', y='fitness', color='srna', hover_data=['condition'], title=title, trendline="ols")
        fig.update_layout(clickmode='event+select')
        fig.update_xaxes(tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1.0], range=[0,1])
        fig.update_traces(marker_size=10)
        return fig
    else:
        return no_data

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
                                html.H3(children="Fitness Finder"),
                            ],
                        ),
                        html.Div(
                            className="container",
                            children=[
                                html.Div(
                                    className="control-block horizontal",
                                    children=[
                                        dcc.Dropdown(
                                            id='dropdown-species',
                                            value='caulo',
                                            clearable=False,
                                            options=[
                                                {'label':'E.Coli', 'value':'ecoli'}, {'label':'C.Crescentus', 'value':'caulo'}
                                            ],
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="control-block horizontal",
                                    children=[
                                        dcc.Dropdown(
                                            id='dropdown-condition',
                                            value=[k for k in description][0],
                                            clearable=False,
                                            placeholder='Select Condition...',
                                            options=[
                                                {'label': value, 'value': key} for key, value in description.items()
                                            ],
                                            multi=True
                                        )
                                    ]
                                ),
                                html.Div(
                                    className="control-block horizontal",
                                    children=[
                                        dcc.Dropdown(
                                            id='dropdown-srna',
                                            value=annotation.iloc[0]['locus_tag'],
                                            clearable=False,
                                            placeholder='Select sRNA...',
                                            options=[
                                                {'label': 'top 5', 'value': 'top'},
                                                {'label': 'worst 5', 'value': 'worst'}
                                                ]
                                        )
                                    ]
                                )
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
                                                dcc.Graph(
                                                    id='graph',
                                                    figure=data_plot([[k for k in description][0]], annotation.iloc[0]['locus_tag'])
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
    Output('graph', 'figure'),
    [Input('dropdown-condition', 'value'),
     Input('dropdown-srna', 'value')])
def update_graph(condition, srna_tag):
    if not isinstance(condition, list): condition = [condition]
    if srna_tag in ['top', 'worst']: return top_ten_plot(condition, srna_tag)
    else: return data_plot(condition, srna_tag)

@app.callback(
    [Output('dropdown-condition', 'options'),
     Output('dropdown-srna', 'value'),
     Output('dropdown-srna', 'options')],
    Input('dropdown-species', 'value'))
def update_options(species):
    global data, description, meta, annotation, groups
    data = all_data[species]['data']
    description = all_data[species]['description']
    meta = all_data[species]['meta']
    annotation = all_data[species]['annotation']
    groups = all_data[species]['groups']
    
    srna_options = [{'label': 'top 5', 'value': 'top'}, {'label': 'worst 5', 'value': 'worst'}] + \
    [{'label': row['locus_tag'], 'value': row['locus_tag']} for i, row in annotation.iterrows() 
     if ((row['start'] <= data['pos']) & (data['pos'] <= row['stop'])).sum() > 0]
    srna_value = 'worst'
    condition_options = [{'label': value, 'value': key} for key, value in description.items()]
    
    return condition_options, srna_value, srna_options
    

def open_browser():
    webbrowser.open_new_tab("http://localhost:8081")

if __name__ == '__main__':
    #open_browser();
    app.run_server(debug=False,port=8081,host='0.0.0.0');