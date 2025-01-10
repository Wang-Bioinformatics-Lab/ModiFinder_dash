import dash
from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_daq as daq
import dash_bootstrap_components as dbc
import pandas as pd
from dash.dash_table.Format import Format, Scheme, Trim
import pickle
import traceback
from tqdm import tqdm
import copy
from furl import furl
import sys
from modifinder.utilities.gnps_types import adduct_mapping

myData = [{'USI1': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906020',
            'USI2': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011905978',
            'Smiles1': 'CCC(C(=O)NCC=CC=C(C)C(C(C)C1C(C(C(O1)C=CC=CC=C(C)C(=O)C2=C(C=CNC2=O)O)O)O)OC)C3(C(C(C(C(O3)C=CC=CC)(C)C)O)O)O',
            'Smiles2': 'CCC(C(=O)NCC=CC=C(C)C(C(C)C1C(C(C(O1)C=CC=CC=C(C)C(=O)C2=C(C=CN(C2=O)C)O)O)O)OC)C3(C(C(C(C(O3)C=CC=CC)(C)C)O)O)O',
            'ppm_tolerance': 40,
            'adduct': '[M+H]1+',
            'name_known':'Kirromycin',
            'name_unknown':'Goldinodox',
            'note':'Example mentioned in the ModiFinder manuscript'}, 
            {'USI1': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906190',
            'USI2': 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906105',
            'Smiles1': 'CC1C=CC=CC=CC(=O)NC2=C(C(=O)C3=C(C2=O)C=C(C(=C3C(=O)C(=CC(C(C(C=CC(CC=C(C(=O)CC1O)C)O)C)O)C)C)O)C)Cl',
            'Smiles2': 'CC1C=CC=CC=C(C(=O)NC2=C(C(=O)C3=C(C2=O)C=C(C(=C3C(=O)C(=CC(C(C(C=CC(CC=C(C(=O)CC1O)C)O)C)O)C)C)O)C)Cl)C',
            'ppm_tolerance': 40,
            'adduct': '[M+H]1+',
            'name_known':'Naphthomycin B',
            'name_unknown':'Naphthomycin A',
            'note':'Example discussed in the ModiFinder manuscript, for peaks 133.0653, 145.0649, and 147.0803, choose the substructures that are result of the fragmentation at the amide bond'},
            {
            'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010118185',
            'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010104042',
            'Smiles1': 'COc1ccc(CC2NCCc3cc(OC)c(OC)cc32)cc1OC',
            'Smiles2': 'COc1ccc(CC2c3cc(OC)c(OC)cc3CCN2C)cc1OC',
            'ppm_tolerance': 40,
            'adduct': '[M+H]1+',
            'name_known':'TETRAHYDROPAPAVERINE',
            'name_unknown':'Laudanosine',
            'note':'An example where ModiFinder performs good.'
            },
            {
            'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010113829',
            'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010125628',
            'Smiles1': 'COC(=O)C1=C(N)Oc2c(c(=O)oc3ccccc23)C1c1ccccc1Br',
            'Smiles2': 'CCOC(=O)C1=C(N)Oc2c(c(=O)oc3ccccc23)C1c1ccccc1Br',
            'ppm_tolerance': 40,
            'adduct': '[M+H]1+',
            'name_known':'methyl 2-amino-4-(2-bromophenyl)-5-oxo-4H-pyrano[3,2-c]chromene-3-carboxylate',
            'name_unknown':'ethyl 2-amino-4-(2-bromophenyl)-5-oxo-4H-pyrano[3,2-c]chromene-3-carboxylate',
            'note':'An example where ModiFinder is not confidant (ambiguity cover issue).'
            },
            {
                'USI1': 'mzspec:GNPS:GNPS-MSMLS:accession:CCMSLIB00005464443',
                'USI2': 'mzspec:GNPS:GNPS-MSMLS:accession:CCMSLIB00005464589',
                'Smiles1': 'N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(O)=O)C(O)=O',
                'Smiles2': 'CCCCCCSCC(C(=O)NCC(=O)O)NC(=O)CCC(C(=O)O)N',
                'ppm_tolerance': 40,
                'adduct': '[M-H]1-',
                'name_known':'GLUTATHIONE REDUCED',
                'name_unknown':'S-HEXYL-GLUTATHIONE',
                'note':'An example for negative adduct'
            },
            {
            'USI1': 'mzspec:GNPS:GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE:accession:CCMSLIB00000847156',
            'USI2': 'mzspec:GNPS:GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE:accession:CCMSLIB00000847119',
            'Smiles1': 'OC1=CC=C(C=C1)C2OC3=C(C(=O)C2C4C(OC5=C(C4=O)C(O)=CC(O)=C5)C6=CC=C(O)C=C6)C(O)=CC(O)=C3',
            'Smiles2': 'COC1=CC=C(C=C1)C2OC3=CC(O)=CC(O)=C3C(=O)C2C4C(OC5=CC(O)=CC(O)=C5C4=O)C6=CC=C(O)C=C6',
            'ppm_tolerance': 40,
            'adduct': '[M+H]1+',
            'name_known':'3-[5,7-dihydroxy-2-(4-hydroxyphenyl)-4-oxo-2,3-dihydrochromen-3-yl]-5,7-dihydroxy-2-(4-hydroxyphenyl)-2,3-dihydrochromen-4-one',
            'name_unknown':'3-[5,7-dihydroxy-2-(4-methoxyphenyl)-4-oxo-2,3-dihydrochromen-3-yl]-5,7-dihydroxy-2-(4-hydroxyphenyl)-2,3-dihydrochromen-4-one',
            'note':'An example where ModiFinder does not perform well. (limited number of shifted peaks)'
            },
            {'USI1': 'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000424840',
            'USI2': 'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000424841',
            'Smiles1': 'CC[C@H](C)[C@H]1C(=O)N2CCC[C@H]2C(=O)O[C@@H](C[C@H](C[C@@H]([C@@H](C3=N[C@H](CS3)/C=C(/C(=O)N[C@H](C(=O)N([C@H](C(=O)N1C)C)C)Cc4ccc(cc4)OC)\C)C)O)C)C(C)(C)C',
            'Smiles2': '',
            'ppm_tolerance': 40,
            'adduct': '[M+H]1+',
            'name_known':'Apratoxin A',
            'name_unknown':'Apratoxin D',
            'note':''},
            ]
column_names = ['USI1', 'USI2', 'Smiles1', 'Smiles2', 'name_known', 'name_unknown', 'note']
hidden_columns = ['USI1', 'USI2', 'Smiles1', 'Smiles2']
supported_adducts =  ['[M+H]1+', '[M-H]1-', '[M+Na]1+', '[M+NH4]1+', '[M+K]1+', '[M+Cl]1-', '[M+Br]1-']
maping = {
    '[M+H]+': '[M+H]1+',
    '[M-H]-': '[M-H]1-',
    '[M+Na]+': '[M+Na]1+',
    '[M+NH4]+': '[M+NH4]1+',
    '[M+K]+': '[M+K]1+',
    '[M+Cl]-': '[M+Cl]1-',
    '[M+Br]-': '[M+Br]1-'
}
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
    for item in ['USI1', 'SMILES1', 'Helpers']:
        inp1.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item + ("' Spectrum IDs (optional)" if item=="Helpers" else ""),id=item, value = args.get(item, ""))],
            style = {'width': '45vw', 'margin': '1vh'}
        ))
    inp2 = [html.H5('Modified Compound',style = {'width': '100%', 'margin': '1vh'})]
    for item in ['USI2', 'SMILES2']:
        inp2.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item + (" (optional)" if item=="SMILES2" else ""),id=item, value = args.get(item, ""))],
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
    for item in ['ppm_tolerance']:#, 'filter_peaks_variable']:
        optionArguments.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
        ))
    
    optionArguments.append(dbc.InputGroup(
            [dbc.InputGroupText('Base Peak Filter Ratio'), dbc.Input(placeholder='filter_peaks_variable',id='filter_peaks_variable', value = args.get('filter_peaks_variable', ""))],
            style={'flex-wrap': 'nowrap', 'minWidth': '300px', 'marginLeft': '1vw'}
        ))
    
    Adduct = args.get('adduct', None)
    if Adduct not in supported_adducts:
        Adduct = adduct_mapping.get(Adduct, None)
        Adduct = maping.get(Adduct, None)


    optionArguments.append(
        dbc.InputGroup(
            [dbc.InputGroupText("Adduct"), dcc.Dropdown(
               supported_adducts, Adduct,
                id = 'adduct', clearable=False, style = {'minWidth': '200px', 'borderTopLeftRadius': 0, 'borderBottomLeftRadius': 0})], style={'flex-wrap': 'nowrap', 'minWidth': '300px', 'marginLeft': '1vw'}
        )
    )

    return html.Div(id = 'inputs', children = 
        [
            # examples
            dbc.Card(
                children=[
                    html.Div(
                        [
                            html.H3("Examples"),
                            dbc.Button([
                                html.I(className="bi bi-chevron-down me-2"), "Open examples"],
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
                dcc.Clipboard(id="input_copy", style={"fontSize":20, "position":"absolute", "right":"10px", "top":"10px"}),

                # error message
                html.Div(id='error-input', style = {'display': 'flex', 'flexDirection': 'row', 'justifyContent': 'center', 'alignItems': 'center', 'width': '100%', 'height': '5vh', 'color': 'red'}),
            ],  style = {'display': 'flex', 'flexDirection': 'column', 'width': '98%', 'padding': '10px', 'justifyContent': 'center', 'alignItems': 'center','margin': '1vw 1vh 1vw 1vh'})
        ], style={'display': 'flex', 'flexDirection': 'column', 'justifyContent': 'center', 'alignItems': 'center', 'width': '100%'})

def get_callbacks(app, helpers):

    @app.callback(
        Output("collapse-button-examples", "children"),
        Output("collapse-examples", "is_open"),
        [Input("collapse-button-examples", "n_clicks"),
         Input("collapse-button-examples", "children")],
        [State("collapse-examples", "is_open")], prevent_initial_call=True
    )
    def toggle_collapse(n, text, is_open):
        button_text = [html.I(className="bi bi-chevron-down me-2"), "Open examples"]
        if "Open examples" in text[1]:
            button_text = [html.I(className="bi bi-chevron-up me-2"), "Close examples"]
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
        State('Helpers', 'value'),
        # State('filter_peaks_method', 'value'),
        State('filter_peaks_variable', 'value'),
        # State('mz_tolerance', 'value'),
        State('ppm_tolerance', 'value'),
        State('adduct', 'value'),
        # State('presence_only', 'on'),
        # State('consider_intensity', 'on'),
        # State('shifted_only', 'on'),
        prevent_initial_call=True
    )
    def _content(n_clicks, data, USI1, USI2, SMILES1, SMILES2, Helpers, filter_peaks_variable, ppm_tolerance, adduct, presence_only = False, consider_intensity=False):
        # remove whitespace from USI
        USI1 = USI1.strip()
        USI2 = USI2.strip()
        try:
            SMILES1 = SMILES1.strip()
        except:
            SMILES1 = None
        
        try:
            SMILES2 = SMILES2.strip()
        except:
            SMILES2 = None

        if SMILES1 == "":
            SMILES1 = None
        
        if SMILES2 == "":
            SMILES2 = None

        temp =  {
            'USI1': USI1,
            'USI2': USI2,
            'SMILES1': SMILES1,
            'SMILES2': SMILES2,
            'Helpers': Helpers,
            # 'filter_peaks_method': filter_peaks_method,
            'filter_peaks_variable': filter_peaks_variable,
            'ppm_tolerance': ppm_tolerance,
            'presence_only': presence_only,
            'consider_intensity': consider_intensity,
            # 'shifted_only': shifted_only,
            'fragmentation_depth':2, 
            "mz_tolerance": 0.01,
            "filter_peaks_method": "intensity",
            # "filter_peaks_variable": 0.01,
            "adduct": adduct
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
        Output('ppm_tolerance', 'value'),
        Output('adduct', 'value')
    ],
    Input('example_table', 'selected_rows'),
    State('example_table', 'data'),
    State('InputData', 'data'), prevent_initial_call=True)
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
                ppm_tolerance = ""
                if 'ppm_tolerance' in data:
                    ppm_tolerance = data['ppm_tolerance']
                Adduct = ""
                if 'adduct' in data:
                    Adduct = data['adduct']
                    if Adduct not in supported_adducts:
                        Adduct = adduct_mapping.get(Adduct, None)
                    Adduct = maping.get(Adduct, None)
                
                return usi1, usi2, smiles1, smiles2, ppm_tolerance, Adduct
            else:
                return "", "", "", "", 40, "[M+H]1+"
        row = local_df[active_rows[0]]
        return row['USI1'], row['USI2'], row['Smiles1'], row['Smiles2'], row['ppm_tolerance'], row['adduct']
    
    @app.callback(
        Output('Helpers', 'value'),
        Input('USI1', 'value'),
        State('USI2', 'value'),
    )
    def update_helpers(usi, usi2):
        specID = None
        if ":" not in usi:
            specID = usi
        elif "accession" in usi:
            specID = usi.split(":")[-1]
        
        if specID is None:
            helper_compounds = []
        else:
            helper_compounds = helpers.get(specID, [])
        
        specID2 = None
        if ":" not in usi2:
            specID2 = usi2
        elif "accession" in usi2:
            specID2 = usi2.split(":")[-1]
        
        if specID2 is not None and specID2 in helper_compounds:
            helper_compounds.remove(specID2)

        return ", ".join(helper_compounds)