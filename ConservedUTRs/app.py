from dash import Dash, no_update
import dash_html_components as html
import dash_core_components as dcc
import dash_table, dash_bio
from dash.dependencies import Input, Output, State

import os
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


# Define Helper Functions

def table_data(df):
    return df[['name', 'length'] + [name for name in df.columns if name.endswith("score")]]

def get_stripped_sequence(sequence):
    return sequence.replace('-', '')

def make_fasta(names, sequences):
    my_fasta = ""
    max_len = 0
    for name, sequence in zip(names, sequences):
        my_fasta += ">{}\n".format(name)
        joined_seq = ''.join(sequence)
        max_len = max(max_len, len(joined_seq))
        my_fasta += "{}\n".format(joined_seq)
    return my_fasta, max_len

def get_inserts(sequence):
    sequence_counter = 0
    insert_counter = 0
    insertions = []
    for nuti in sequence:
        if nuti != '-': 
            if insert_counter != 0:
                insertions.append((sequence_counter, insert_counter))
                insert_counter = 0
            sequence_counter += 1
        else:
            insert_counter += 1
    if insert_counter != 0: insertions.append((sequence_counter, insert_counter))
    return insertions

def get_sequence_with_inserts(sequence, inserts, absolute_positions=False):
    my_seq = ""
    if not inserts: return sequence
    my_inserts = inserts[::-1]
    current_insert = my_inserts.pop()
    current_sequence_position = 0
    slice_position = 0
    for i,n in enumerate(sequence):
        if current_sequence_position == current_insert[0]:
            my_seq += sequence[slice_position:i] + '-' * current_insert[1]
            slice_position = i
            if my_inserts: current_insert = my_inserts.pop()
            else: return my_seq + sequence[i:]
        if absolute_positions: current_sequence_position += 1
        elif n != '-': current_sequence_position += 1
    my_seq += sequence[slice_position:]
    if current_insert[0] == current_sequence_position: my_seq += '-' * current_insert[1]
    return my_seq

def get_absolute_inserts(sequence, inserts):
    my_inserts = np.array(inserts)
    seq_ins = get_inserts(sequence)
    for seq_in in seq_ins: my_inserts[my_inserts[:,0] >= seq_in[0]] += seq_in[1]
    return [tuple(seq_in) for seq_in in my_inserts]
    
        
def multialigned(seqname, alignments):
    if not alignments: return make_fasta(['no_data'], [''])
    if len(alignments) == 1: return make_fasta([seqname, alignments[0][0]], alignments[0][1:3])
    spec, seqaln, refaln, score = alignments[0]
    species = [spec]
    refalns = [refaln]
    for current_species, current_seqaln, current_refaln, _ in alignments[1:]:
        seq_inserts = get_inserts(seqaln)
        current_seq_inserts = get_inserts(current_seqaln)
        abs_current_seqaln_inserts = get_absolute_inserts(current_seqaln, seq_inserts)
        abs_new_seqaln_inserts = get_absolute_inserts(seqaln, current_seq_inserts)
        new_refaln = get_sequence_with_inserts(refaln, abs_current_seqaln_inserts, absolute_positions=True)
        seqaln = get_sequence_with_inserts(seqaln, seq_inserts)
        for refaln in refalns: refaln = get_sequence_with_inserts(refaln, abs_new_seqaln_inserts, absolute_positions=True)
        refalns.append(new_refaln)
        species.append(current_species)
    return make_fasta([seqname] + species, [seqaln] + refalns)
    
alignments = [
    ['test1', 'ABCDE-FG', 'ABCDEFFG', 1],
    ['test2', 'ABCDEFG', 'ABC-EFG', 1],
    ['test1', 'ABC-DEFG', 'ABCCDEFG', 1],
    ]
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
                                                    sort_action="native",
                                                    filter_action='native',
                                                    page_size=100,
                                                    column_selectable='multi',
                                                    style_header=
                                                    {
                                                    'fontWeight': 'bold',
                                                    'border': 'thin lightgrey solid',
                                                    'backgroundColor': 'rgb(100, 100, 100)',
                                                    'color': 'white'
                                                    },
                                                    style_cell={
                                                    'fontFamily': 'Open Sans',
                                                    'textAlign': 'left',
                                                    'width': '80px',
                                                    'minWidth': '80px',
                                                    'maxWidth': '80px',
                                                    'whiteSpace': 'no-wrap',
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis',
                                                    'backgroundColor': 'Rgb(230,230,250)'
                                                    },
                                                    style_data_conditional=[
                                                    {
                                                    'if': {'row_index': 'odd'},
                                                    'backgroundColor': 'rgb(248, 248, 248)'
                                                    },
                                                    {
                                                    'if': {'column_id': 'name'},
                                                    'color': 'black',
                                                    'fontWeight': 'bold',
                                                    'textAlign': 'center'
                                                    }
                                                    ],
                                                    fixed_rows={'headers': True, 'data': 0}
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
                                                    colorscale={'A':'#11BB11', 'T':'#CC0000', 'C':'#1111CC', 'G':'#FF8011', '-': '#ffffff'},
                                                    textsize=12,
                                                    tilewidth=30,
                                                    overview='none',
                                                    height='180px',
                                                    showconsensus=False
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
                                                                {'label': i[:-6], 'value': i[:-6]} for i in table_data(df).columns if (i.endswith("score") and not i.startswith("average"))
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

@app.callback([Output('my-alignment-viewer', "data"),
               Output('my-alignment-viewer', 'numtiles'),
               Output('my-alignment-viewer', 'showconsensus')],
              [Input('table', "active_cell"),
               Input('species-checklist', "value")],
              state=State('table', "page_current"))
def get_active_alignment(active_cell, checked_species, page_current):
    if active_cell:
        if page_current is None: add = 0
        else: add = page_current * 100
        row_index = add + active_cell['row']
        if active_cell['column_id'] not in ['name', 'length', 'average_score']:
            my_fasta = ">{}\n{}-\n".format(df.iloc[row_index]['name'], df.iloc[row_index][active_cell['column_id'][0:-6] + "_seqaln"])
            my_fasta += ">{}\n{}-".format(active_cell['column_id'][0:-6], df.iloc[row_index][active_cell['column_id'][0:-6] + "_refaln"])
            length = len(df.iloc[row_index][active_cell['column_id'][0:-6] + "_refaln"])
            consens = False
        elif active_cell['column_id'] in ['name', 'average_score']:
            if checked_species is None: return no_update, no_update, no_update
            seqname = df.iloc[row_index]['name']
            alignments = sorted([[species] + list(df[[species + "_seqaln", species + "_refaln", species + "_score"]].iloc[row_index]) for species in checked_species], key=lambda x:x[3], reverse=True)
            my_fasta, length = multialigned(seqname, alignments)
            #sequences = [df.iloc[row_index]['sequence']] + [df.iloc[row_index][species + "_refaln"] for species in checked_species]
            #aligned = mult_align(sequences, gop=-100, scale=3.0)
            #names = [seqname] + [species for species in checked_species]
            #my_fasta, length = make_fasta(names, sequences)
            consens = True
        else:
            return no_update, no_update, no_update
        return my_fasta, length, consens
    return no_update, no_update, no_update

if __name__ == '__main__':
    app.run_server(debug=True,port=8083,host='0.0.0.0');
