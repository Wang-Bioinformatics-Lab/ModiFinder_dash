from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_daq as daq
import dash_bootstrap_components as dbc
import pandas as pd
from dash.dash_table.Format import Format, Scheme, Trim
import pickle
import traceback
from tqdm import tqdm
import copy

myData = [{'USI1': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906020',
            'USI2': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011905978',
            'Smiles1': 'CCC(C(=O)NCC=CC=C(C)C(C(C)C1C(C(C(O1)C=CC=CC=C(C)C(=O)C2=C(C=CNC2=O)O)O)O)OC)C3(C(C(C(C(O3)C=CC=CC)(C)C)O)O)O',
            'Smiles2': 'CCC(C(=O)NCC=CC=C(C)C(C(C)C1C(C(C(O1)C=CC=CC=C(C)C(=O)C2=C(C=CN(C2=O)C)O)O)O)OC)C3(C(C(C(C(O3)C=CC=CC)(C)C)O)O)O',
            'name_known':'Kirromycin',
            'name_unknown':'Goldinodox',
            'note':'Example discussed in ModiFinder paper'}, 
            {'USI1': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906190',
            'USI2': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906105',
            'Smiles1': 'CC1C=CC=CC=CC(=O)NC2=C(C(=O)C3=C(C2=O)C=C(C(=C3C(=O)C(=CC(C(C(C=CC(CC=C(C(=O)CC1O)C)O)C)O)C)C)O)C)Cl',
            'Smiles2': 'CC1C=CC=CC=C(C(=O)NC2=C(C(=O)C3=C(C2=O)C=C(C(=C3C(=O)C(=CC(C(C(C=CC(CC=C(C(=O)CC1O)C)O)C)O)C)C)O)C)Cl)C',
            'name_known':'Naphthomycin B',
            'name_unknown':'Naphthomycin A',
            'note':'Example discussed in ModiFinder paper'},
            ]
column_names = ['USI1', 'USI2', 'Smiles1', 'Smiles2', 'name_known', 'name_unknown', 'note']
hidden_columns = ['USI1', 'USI2', 'Smiles1', 'Smiles2']
columns = [{"name": i, "id": i} for i in column_names]

examples = dash_table.DataTable(
        id='example_table',
        columns=columns,
        data=myData,
        style_cell={
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'textAlign': 'left',
        'maxWidth': '25%',
        },
        style_data = {        'whiteSpace': 'normal',
        'height': 'auto',},
        tooltip_delay=0,
        tooltip_duration=None,
        tooltip_header={i['name']: i['name'] for i in columns},
        page_current=0,
        row_selectable="single",
        style_table={
                # 'maxHeight': '50ex',
                'overflowY': 'scroll',
                'width': '100%',
                'minWidth': '100%',
            },

        style_cell_conditional=[
            {
            'if': {'column_id': c},
            'display': 'none'
            } for c in hidden_columns
        ]
    )

def get_layout(passedArgs):
    
    args = copy.deepcopy(passedArgs)
        
    inp1 = [html.H5('Known Compound',style = {'width': '100%', 'margin': '1vh'})]
    for item in ['USI1', 'SMILES1']:
        inp1.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
            style = {'width': '45vw', 'margin': '1vh'}
        ))
    inp2 = [html.H5('Modified Compound',style = {'width': '100%', 'margin': '1vh'})]
    for item in ['USI2', 'SMILES2']:
        inp2.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
            style = {'width': '45vw', 'margin': '1vh'}
        ))
    layer1 = []
    layer1.append(html.Div(children = inp1, style = {'border': '1px solid #aabbdd', 'display': 'flex', 'flexDirection': 'column', 'justifyContent': 'center', 'alignItems': 'center', 'padding':'10px'}))
    layer1.append(html.Div(children = inp2, style = {'border': '1px solid #aabbdd', 'display': 'flex', 'flexDirection': 'column', 'justifyContent': 'center', 'alignItems': 'center', 'padding':'10px'}))
    layer1 = html.Div(children = layer1, style = {'display': 'flex', 'flexDirection': 'row', 'flexWrap': 'wrap', 'justifyContent': 'space-around', 'width': '100%', 'minHeight': '5vh'})

    myOptions = []
    options = []#[{'label': 'Presence Only', 'value': 'presence_only'}, {'label': 'Consider Intensity', 'value': 'consider_intensity'}]#, {'label': 'Shifted Only', 'value': 'shifted_only'}]
    for item in options:
        myOptions.append(daq.BooleanSwitch(
            id=item['value'],
            on=args.get(item['value'], False),
            label=item['label'],
            labelPosition="Left" , style={'display': 'inline-block', 'margin': '2vh 1vw'})
        )
    
    switches = None
    if len(myOptions) > 0:
        switches = html.Div(id='options',
                        children=myOptions,
                    style={'display': 'flex', 'flexDirection': 'row','margin-right': '1vw'}),

    optionArguments = []
    # optionArguments = [dcc.Dropdown(['intensity', 'top_k', 'none'], 'top_k', id='filter_peaks_method')]
    for item in ['ppm']:#, 'filter_peaks_variable']:
        optionArguments.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
        ))

    return html.Div(id = 'inputs', children = 
        [
            # examples
            dbc.Card(
                children=[
                    html.Div(
                        [
                            html.H3("Examples"),
                            dbc.Button(
                                "Open examples",
                                id="collapse-button-examples",
                                className="mb-3",
                                color="primary",
                                n_clicks=0,
                            ),
                        ], style={'display': 'flex', 'flexDirection': 'row', 'justifyContent': 'space-between', 'alignItems': 'start'}
                    ),
                    dbc.Collapse(
                        examples,
                        id="collapse-examples",
                        is_open=False,
                    ),
                ], style={"width": "98%", 'margin': '1vw 1vh 1vw 1vh', 'padding': '10px'}
            ),

            dbc.Card(children=[
                html.H3("Input", style={'alignSelf': 'flex-start'}),
                
                # compound inputs
                layer1,
                
                # options
                html.Div(children = [
                        switches,
                        html.Div(id = 'arguments', children = optionArguments, style = {'display': 'flex', 'flexDirection': 'row', 'justifyContent': 'center', 'minHeight': '5vh', 'alignItems': 'center'}),
                    ], style = {'display': 'flex', 'flexDirection': 'row', 'justifyContent': 'center', 'alignItems': 'center'}),
                
                
                # update button
                dbc.Button('update', id='update', style = {'width': '300px', 'margin': 'auto'}),
            ],  style = {'display': 'flex', 'flexDirection': 'column', 'width': '98%', 'padding': '10px', 'justifyContent': 'center', 'alignItems': 'center','margin': '1vw 1vh 1vw 1vh'})
        ], style={'display': 'flex', 'flexDirection': 'column', 'justifyContent': 'center', 'alignItems': 'center', 'width': '100%'})

def get_callbacks(app):
    @app.callback(
        Output("collapse-button-examples", "children"),
        Output("collapse-examples", "is_open"),
        [Input("collapse-button-examples", "n_clicks"),
         Input("collapse-button-examples", "children")],
        [State("collapse-examples", "is_open")], prevent_initial_call=True
    )
    def toggle_collapse(n, text, is_open):
        button_text = "Open examples"
        if text == "Open examples":
            button_text = "Close examples"
        if n:
            return button_text, not is_open
        return button_text, is_open

    @app.callback(
        Output('inputs', 'children'),
        Input('urlData', 'data'),
        prevent_initial_call=True
    )
    def _content(data):
        return get_layout(data)

    @app.callback(
        Output('InputData', 'data'),
        Input('update', 'n_clicks'),
        State('InputData', 'data'),
        State('USI1', 'value'),
        State('USI2', 'value'),
        State('SMILES1', 'value'),
        State('SMILES2', 'value'),
        # State('filter_peaks_method', 'value'),
        # State('filter_peaks_variable', 'value'),
        # State('mz_tolerance', 'value'),
        State('ppm', 'value'),
        # State('presence_only', 'on'),
        # State('consider_intensity', 'on'),
        # State('shifted_only', 'on'),
        prevent_initial_call=True
    )
    def _content(n_clicks, data, USI1, USI2, SMILES1, SMILES2, ppm, presence_only = False, consider_intensity=False):
        # remove whitespace from USI
        USI1 = USI1.strip()
        USI2 = USI2.strip()
        SMILES1 = SMILES1.strip()
        SMILES2 = SMILES2.strip()
        temp =  {
            'USI1': USI1,
            'USI2': USI2,
            'SMILES1': SMILES1,
            'SMILES2': SMILES2,
            'args':{
                # 'filter_peaks_method': filter_peaks_method,
                # 'filter_peaks_variable': filter_peaks_variable,
                'ppm': ppm,
                'presence_only': presence_only,
                'consider_intensity': consider_intensity,
                # 'shifted_only': shifted_only,
                'fragmentation_depth':3, 
                "mz_tolerance": 0.1,
                "filter_peaks_method": "top_k",
                "filter_peaks_variable": 30
            }
        }
        if data is None:
            data = {}
        data.update(temp)
        return data
    @app.callback(
    [
        Output('USI1', 'value'),
        Output('USI2', 'value'),
        Output('SMILES1', 'value'),
        Output('SMILES2', 'value'),
    ],
    Input('example_table', 'selected_rows'),
    State('example_table', 'data'),
    State('InputData', 'data'))
    def update_inputs_from_table(active_rows, local_df, data):
        if active_rows is None or len(active_rows) == 0:
            if data is not None:
                usi1 = ""
                if 'USI1' in data:
                    usi1 = data['USI1']
                usi2 = ""
                if 'USI2' in data:
                    usi2 = data['USI2']
                smiles1 = ""
                if 'SMILES1' in data:
                    smiles1 = data['SMILES1']
                smiles2 = ""
                if 'SMILES2' in data:
                    smiles2 = data['SMILES2']
                return usi1, usi2, smiles1, smiles2
            else:
                return "", "", "", ""
        row = local_df[active_rows[0]]
        return row['USI1'], row['USI2'], row['Smiles1'], row['Smiles2']