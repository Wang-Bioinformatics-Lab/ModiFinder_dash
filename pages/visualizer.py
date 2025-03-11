import dash
from dash import Dash, html, dcc, Input, Output, State, callback, ctx
import dash_bootstrap_components as dbc
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
import numpy as np
import base64
from io import BytesIO
import sys
from PIL import Image
import urllib
# sys.path.append('../')
# sys.path.append('../SmallMol_Mod_Site_Localization')
# from SmallMol_Mod_Site_Localization.modifinder.utilities import visualizer as vis
# from SmallMol_Mod_Site_Localization.modifinder import alignment as align
# from SmallMol_Mod_Site_Localization.modifinder.utilities.network import get_spectrum
# import SmallMol_Mod_Site_Localization.modifinder.utilities.spectra_utils as su
# import SmallMol_Mod_Site_Localization.modifinder.utilities.mol_utils as mu

import modifinder.utilities.visualizer as mf_vis
import modifinder.utilities.mol_utils as mu
from modifinder import Spectrum
from modifinder.engines.alignment.CosineAlignmentEngine import _cosine_fast

from flask import Flask, send_file, request, jsonify
import json
from app import app
from furl import furl
from myopic_mces import MCES

# dash.register_page(__name__)
# app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP])
dash_app = Dash(
    name="visualizer dashboard",
    server=app,
    url_base_pathname="/visualizer/",
    external_stylesheets=[dbc.themes.BOOTSTRAP],
)

dash_app.index_string = """<!DOCTYPE html>
<html>
    <head>
        <!-- Umami Analytics -->
        <script async defer data-website-id="88fbdd25-49e4-4b21-be94-689631087d04" src="https://analytics.gnps2.org/umami.js"></script>
        <script async defer data-website-id="74bc9983-13c4-4da0-89ae-b78209c13aaf" src="https://analytics.gnps2.org/umami.js"></script>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>"""

button_tab_style = {
    'background-color': 'rgba(0,0,0,0.1)',
    'color': 'rgba(0,0,0,0.6)',
    'cursor': 'pointer',
    'font-size': '1.5rem',
    'margin': '10px',
    'flex': '1',
    'margin': 0,
    'border-left': 'none',
    'border-right': 'none',
    'border-top': 'none',
    'border-bottom': '1px solid blue',
}

button_tab_active_style = {
    'background-color': 'transparent',
    'color': 'black',
    'cursor': 'pointer',
    'font-size': '1.5rem',
    'margin': '10px',
    'flex': '1',
    'border-left': '1px solid blue',
    'border-right': '1px solid blue',
    'border-top': '1px solid blue',
    'border-bottom': 'none',
    'margin': 0,
}


NAVBAR = dbc.Navbar(
    children=[
        dbc.NavbarBrand(
            html.Img(src="https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/img/logo/GNPS2_logo_blue-grey-black.png", width="120px"),
            href="https://www.cs.ucr.edu/~mingxunw/"
        ),
        dbc.Nav(
            [
                dbc.NavItem(dbc.NavLink("Wang Bioinformatics Lab", href="https://www.cs.ucr.edu/~mingxunw/")),
                dbc.NavItem(dbc.NavLink("ModiFinder Dashboard", href="https://modifinder.gnps2.org/")),
                dbc.NavItem(dbc.NavLink("Visualizer", href="https://modifinder.gnps2.org/visualizer")),
                dbc.NavItem(dbc.NavLink("Tutorial", href="https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/modifinder-web/")),
                dbc.NavItem(dbc.NavLink("Documentation", href="https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/modifinder/")),
            ],
        navbar=True)
    ],
    color="light",
    dark=False,
    sticky="top",
    style={"padding": "10px 20px"},
)


DRAW_MOLECULE =  html.Div([
    html.H1("Molecule Drawer"),
    dbc.Card([
    dbc.InputGroup(
            [dbc.InputGroupText('Smiles1'), dbc.Input(placeholder='SMILES, InChI, Spectrum ID, or USI',id='Mol1', value = "")],
            style = {'width': '90vw', 'margin': '1vh auto'}
    ),
    dbc.InputGroup(
            [dbc.InputGroupText('Smiles2'), dbc.Input(placeholder='SMILES, InChI, Spectrum ID, or USI',id='Mol2', value = "")],
            style = {'width': '90vw', 'margin': '1vh auto'}
    ),
    dbc.Checklist(
        id='molecule-boolean-inputs',
        options=[
            {'label': 'Highlight common', 'value': 'highlight_common'},
            {'label': 'Highlight removed', 'value': 'highlight_removed'},
            {'label': 'Highlight added', 'value': 'highlight_added'},
            {'label': 'Modification only', 'value': 'modification_only'},
            {'label': 'Show Legend', 'value': 'show_legend'},
        ],
        value=['highlight_common', 'highlight_removed', 'highlight_added', 'show_legend'],
        switch=True,
        style={'display': 'flex', 'flex-direction': 'row', 'justify-content': 'space-around', 'margin': '1vh auto', 'width': '80vw', 'wrap': 'wrap'}
    ),
    dcc.Clipboard(id='draw_molecule_clipboard', style={'position': 'absolute', 'top': '5px', 'right': '5px'}),]),
    dcc.Loading(html.Div(id='output-molecule')),

    ], style={'margin': 'auto', 'width': '95vw'})

DRAW_ALIGNMENT = html.Div([
    html.H1("Alignment Drawer"),
    dbc.Card([
    dbc.InputGroup(
            [dbc.InputGroupText('Spec1'), dbc.Input(placeholder='Spectrum ID or USI',id='Spec1', value = "")],
            style = {'width': '90vw', 'margin': '1vh auto'}
    ),
    dbc.InputGroup(
            [dbc.InputGroupText('Spec2'), dbc.Input(placeholder='Spectrum ID or USI',id='Spec2', value = "")],
            style = {'width': '90vw', 'margin': '1vh auto'}
    ),
    dbc.Checklist(
        id='spectra-boolean-inputs',
        options=[
            {'label': 'Normalize Peaks', 'value': 'normalize_peaks'},
            {'label': 'Draw Mapping Lines', 'value': 'draw_mapping_lines'},
            {'label': 'Flipped Representation', 'value': 'flipped'},
        ],
        value=['normalize_peaks', 'draw_mapping_lines', 'flipped'],
        switch=True,
        style={'display': 'flex', 'flex-direction': 'row', 'justify-content': 'space-around', 'margin': '1vh auto', 'width': '80vw', 'wrap': 'wrap'}
    ),
    dcc.Clipboard(id='draw_spectra_clipboard', style={'position': 'absolute', 'top': '5px', 'right': '5px'}),]),
    dcc.Loading(html.Div(id='output-spectra'))
    ], style={'margin': 'auto', 'width': '95vw'})


CONTRIBUTORS_DASHBOARD = [
    dbc.CardHeader(html.H5("Contributors")),
    dbc.CardBody(
        [
            "Reza Shahneh - UC Riverside",
            html.Br(),
            html.Br(),
            html.H5("Citation"),
            html.A('Mohammad Reza Zare Shahneh, Michael Strobel, Giovanni Andrea Vitale, Christian Geibel, Yasin El Abiead, Neha Garg, Berenike Wagner, Karl Forchhammer, Allegra Aron, Vanessa V Phelan, Daniel Petras, Mingxun Wang. "ModiFinder: Tandem Mass Spectral Alignment Enables Structural Modification Site Localization" Journal of the American Society for Mass Spectrometry, doi:10.1021/jasms.4c00061, PMID: 38830143', 
                    href="https://doi.org/10.1021/jasms.4c00061"),
            html.Br(),
            html.Br(),
            html.A('Checkout our other work!', 
                href="https://www.cs.ucr.edu/~mingxunw/")
        ]
    )
]

modification_options = ['highlight_common', 'highlight_removed', 'highlight_added', 'modification_only', 'show_legend']
spectra_options = ['normalize_peaks', 'draw_mapping_lines', 'flipped']

dash_app.layout = html.Div([
    NAVBAR,
    html.Div(children=[
        html.Div(id='tabs', children=[
            html.Button('Molecule', id='molecule-tab-button', style={'display': 'none'}),
            html.Button('Spectra', id='spectra-tab-button', style={'display': 'none'}),
        ], style={'display': 'flex', 'flex-direction': 'row', 'width': '100%', 'justify-content': 'center', 'align-items': 'center'}),
        html.Div(id='molecule-tab', children=DRAW_MOLECULE, style={'display': 'none'}),
        html.Div(id='spectra-tab', children=DRAW_ALIGNMENT, style={'display': 'none'}),
    ]),
    dbc.Card(CONTRIBUTORS_DASHBOARD, style={"width": "98%", "margin": "auto", 'align-self': 'flex-end', 'margin-bottom': '10px'}),
    dcc.Store(id='selected_tab', data='Molecule'),
    dcc.Location(id='url', refresh=False),
], style={'display': 'flex', 'flex-direction': 'column', 'min-height': '100vh'})

def convert_stats(stats_array):
    stats_section = []
    for stat in stats_array:
        stats_section.append(
            html.Div(
                children=[
                    dbc.Badge(stat[0], color="info", className="me-1"),
                    html.P(
                        children=stat[1],
                        style={"margin": "0px", "margin-left": "5px"},
                    ),
                ],
                style={
                    "display": "flex",
                    "flexDirection": "row",
                    "alignItems": "center",
                    "justifyContent": "space-between",
                    "backgroundColor": "#f8f9fa",
                    "margin": "3px",
                    "padding": "0px",
                    "borderRadius": "5px",
                },
            )
        )
    return stats_section

def png_to_showable_src(png):
    if png.max() <= 1:
        png = (png*255).astype(np.uint8)
    img = Image.fromarray(png)
    buffer = BytesIO()
    img.save(buffer, format='PNG')
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    return f'data:image/png;base64,{img_base64}'

def get_stats(molecule, draw_mol_png = None):
    h_added = Chem.AddHs(molecule)
    stats = [
        ["Molecular Weight", rdMolDescriptors.CalcExactMolWt(molecule)],
        ["Molecular Formula", rdMolDescriptors.CalcMolFormula(molecule)],
        ["Number of Atoms", h_added.GetNumAtoms()],
        ["Number of Heavy Atoms", h_added.GetNumHeavyAtoms()],
        ["Number of Bonds", molecule.GetNumBonds()],
        ["Number of Rings", molecule.GetRingInfo().NumRings()],
        ['SMILES', Chem.MolToSmiles(molecule)],
        ['InChI', Chem.MolToInchi(molecule)],
        ['InChI Key', Chem.InchiToInchiKey(Chem.MolToInchi(molecule))],
    ]
    stats_section = []

    if draw_mol_png is not None:
        if type(draw_mol_png) == str:
            img = draw_mol_png
        else:
            img = png_to_showable_src(draw_mol_png)
        stats_section.append(html.Img(src=img, style={'margin': 'auto', 'height': '20vh'}))

    stats_section += convert_stats(stats)
    stats = dbc.Card(
        children=stats_section,
        style={"flex": "1", "width": "40%", "margin": "10px"},
    )
    return stats


@dash_app.callback(
    Output('selected_tab', 'data'),
    Input('molecule-tab-button', 'n_clicks'),
    Input('spectra-tab-button', 'n_clicks'),
    Input('url', 'href'),
    State('selected_tab', 'data'))
def update_selected_tab(molecule_clicks, spectra_clicks, href, selected_tab):
    ctx = dash.callback_context
    if not ctx.triggered:
        return selected_tab
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    f = furl(href)
    arguments = f.args
    url_tab = arguments.get('selected_tab', selected_tab)
    if button_id == 'molecule-tab-button':
        return 'Molecule'
    elif button_id == 'spectra-tab-button':
        return 'Spectra'
    else:
        return url_tab
    

@dash_app.callback(
    Output('molecule-tab', 'style'),
    Output('spectra-tab', 'style'),
    Output('molecule-tab-button', 'style'),
    Output('spectra-tab-button', 'style'),
    Input('selected_tab', 'data')
)
def update_tabs(selected_tab):
    if selected_tab == 'Molecule':
        return {'display': 'block'}, {'display': 'none'}, button_tab_active_style, button_tab_style
    else:
        return {'display': 'none'}, {'display': 'block'}, button_tab_style, button_tab_active_style

def draw_molecule(Mol1, Mol2, **kwargs):
    png = mf_vis.draw_modifications(Mol1, Mol2, **kwargs)
    # png is numpy array, convert to base64
    img = Image.fromarray(png)
    buffer = BytesIO()
    img.save(buffer, format='PNG')
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    return f'data:image/png;base64,{img_base64}'

@dash_app.callback(
    Output('output-molecule', 'children'),
    Input('Mol1', 'value'),
    Input('Mol2', 'value'),
    Input('molecule-boolean-inputs', 'value'),
)
def draw_molecule_output(Mol1, Mol2, boolean_inputs):
    try:
        if Mol1 == '':
            Mol1 = None
        if Mol2 == '':
            Mol2 = None
        stats = None
        if Mol1 is None or Mol2 is None:
            if Mol1 is None and Mol2 is None:
                return None
            if Mol1 is None:
                input = Mol2
            else:
                input = Mol1
            input_mol = mu._get_molecule(input)
            png = mf_vis.draw_molecule(input_mol)
            stats = get_stats(input_mol, draw_mol_png = png_to_showable_src(png))
            return stats
        else:
            kwargs = {k: True for k in boolean_inputs}
            for option in modification_options:
                if option not in kwargs:
                    kwargs[option] = False
            png = mf_vis.draw_modifications(Mol1, Mol2, **kwargs)
            
            mol1, mol2 = mu._get_molecules(Mol1, Mol2)
            smiles1 = Chem.MolToSmiles(mol1)
            smiles2 = Chem.MolToSmiles(mol2)
            mces_value = MCES(smiles1, smiles2)
            print(mces_value)
            result = mu.get_transition(mol1, mol2)
            added = result['modified_added_edges_inside'] + result['modified_added_edges_bridge']
            removed = result['modified_removed_edges_inside'] + result['modified_removed_edges_bridge']
            fpgen = AllChem.GetMorganGenerator(radius=2)
            fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in [mol1, mol2]]
            tanimoto = round(DataStructs.TanimotoSimilarity(fps[0],fps[1]), 2)
            edit_stats = [
                ["Edit Distance", len(added) + len(removed)],
                ["Number of Added", len(added)],
                ["Number of Removed", len(removed)],
                ["Tanimoto Similarity (Morgan, radius=2)", tanimoto],
                ["myopic MCES distance", round(mces_value[1], 2)],
            ]

            edit_stats_card = dbc.Card(
                children=convert_stats(edit_stats),
                style={"width": "20%", "padding": "10px"},
            )

            copy_mol1 = Chem.Mol(mol1)
            copy_mol2 = Chem.Mol(mol2)
            mcs1 = rdFMCS.FindMCS([copy_mol1, copy_mol2])
            mcs_mol = Chem.MolFromSmarts(mcs1.smartsString)
            mol1_image_png = png_to_showable_src(mf_vis.draw_modifications(mol1, mcs_mol, **kwargs))
            mol2_image_ong = png_to_showable_src(mf_vis.draw_modifications(mcs_mol, mol2, **kwargs))

            stats = html.Div([get_stats(mu._get_molecule(Mol1), draw_mol_png = mol1_image_png), get_stats(mu._get_molecule(Mol2), draw_mol_png = mol2_image_ong)],
                                style={'display': 'flex', 'flex-direction': 'row', 'justify-content': 'space-around', 'align-items': 'center'})

            img_base64 = png_to_showable_src(png)
            return [html.Div([html.Div(["Merged Image", html.Img(src=img_base64)], 
                                       style={'display': 'flex', 'flex-direction': 'column', 'justify-content': 'center', 'align-items': 'center', 'flex':'1'}),
                               edit_stats_card],
                            style={'display': 'flex', 'flex-direction': 'row', 'justify-content': 'space-around', 'align-items': 'center'}),
                                stats]
    except Exception as e:
        return str(e)

@dash_app.callback(
    Output('output-spectra', 'children'),
    Input('Spec1', 'value'),
    Input('Spec2', 'value'),
    Input('spectra-boolean-inputs', 'value'),
)
def update_spectra_output(Spec1, Spec2, boolean_inputs):
    try:
        if Spec1 == '':
            Spec1 = None
        if Spec2 == '':
            Spec2 = None

        kwargs = {k: True for k in boolean_inputs}
        for option in spectra_options:
            if option not in kwargs:
                kwargs[option] = False
        if Spec1 is None or Spec2 is None:
            if Spec1 is None and Spec2 is None:
                return None
            if Spec1 is None:
                input = Spec2
            else:
                input = Spec1
            input = Spectrum(input, ignore_adduct_format=True)
            png = mf_vis.draw_spectrum(input, **kwargs)
        else:
            Spec1 = Spectrum(Spec1, ignore_adduct_format=True)
            Spec2 = Spectrum(Spec2, ignore_adduct_format=True)
            cosine, matches = _cosine_fast(Spec1, Spec2, 0.1, 40, True)
            png = mf_vis.draw_alignment([Spec1, Spec2], [matches], **kwargs)
        
        img = png_to_showable_src(png)
        return html.Img(src=img, style={'margin': 'auto', 'height': '50vh'})
    except Exception as e:
        return str(e)


@dash_app.callback(
    Output('draw_molecule_clipboard', 'content'),
    Input('draw_molecule_clipboard', 'n_clicks'),
    State('Mol1', 'value'),
    State('Mol2', 'value'),
    State('molecule-boolean-inputs', 'value'),
    State('url', 'href'),
    State('selected_tab', 'data'))
def update_clipboard_molecule(n_clicks, Mol1, Mol2, boolean_inputs, href, selected_tab):
    f = furl(href)
    arguments = {}
    arguments['Mol1'] = Mol1
    arguments['Mol2'] = Mol2
    arguments.update({k: True for k in boolean_inputs})
    for option in modification_options:
        if option not in arguments:
            arguments[option] = False
    arguments['selected_tab'] = selected_tab
    f.args = arguments
    return f.url

@dash_app.callback(
    Output('draw_spectra_clipboard', 'content'),
    Input('draw_spectra_clipboard', 'n_clicks'),
    State('Spec1', 'value'),
    State('Spec2', 'value'),
    State('spectra-boolean-inputs', 'value'),
    State('url', 'href'),
    State('selected_tab', 'data'))
def update_clipboard_spectra(n_clicks, Spec1, Spec2, boolean_inputs, href, selected_tab):
    f = furl(href)
    arguments = {}
    arguments['Spec1'] = Spec1
    arguments['Spec2'] = Spec2
    arguments.update({k: True for k in boolean_inputs})
    for option in spectra_options:
        if option not in arguments:
            arguments[option] = False
    arguments['selected_tab'] = selected_tab
    f.args = arguments
    return f.url


# fill the input fields with the URL arguments
@dash_app.callback(
    Output('Mol1', 'value'),
    Output('Mol2', 'value'),
    Output('molecule-boolean-inputs', 'value'),
    Output('Spec1', 'value'),
    Output('Spec2', 'value'),
    Output('spectra-boolean-inputs', 'value'),
    Input('url', 'href'),
    State('molecule-boolean-inputs', 'value'),
    State('spectra-boolean-inputs', 'value'),
    allow_duplicate=True
)
def update_input_fields(href, modification_options, spectra_options):
    f = furl(href)
    arguments = f.args
    Mol1 = arguments.get('Mol1', '')
    Mol2 = arguments.get('Mol2', '')
    Spec1 = arguments.get('Spec1', '')
    Spec2 = arguments.get('Spec2', '')
    molecule_boolean_inputs = []
    spectra_boolean_inputs = []
    for key in arguments:
        if key in modification_options and arguments[key] == 'True':
            molecule_boolean_inputs.append(key)
        elif key in spectra_options and arguments[key] == 'True':
            spectra_boolean_inputs.append(key)
    if len(molecule_boolean_inputs) == 0:
        molecule_boolean_inputs = modification_options
    if len(spectra_boolean_inputs) == 0:
        spectra_boolean_inputs = spectra_options
    return Mol1, Mol2, molecule_boolean_inputs, Spec1, Spec2, spectra_boolean_inputs


@app.route("/api/visualizer/<function_name>", methods=['GET'])
def api(function_name):
    vis_functions = mf_vis.return_public_functions()
    if function_name not in vis_functions:
        return jsonify({'error': 'Function not found'}), 404
    try:
        args = request.args
        kwargs = {"ignore_adduct_format": True}
        file_type = 'png'
        print(args)
        for key in args:
            if '.png' in args[key]:
                kwargs[key] = args[key].split('.png')[0]
                file_type = 'png'
            elif '.svg' in args[key]:
                kwargs[key] = args[key].split('.svg')[0]
                file_type = 'svg'
            else:
                kwargs[key] = args[key]
            if key == 'output_type':
                file_type = args[key]
        
        # all args are strings, convert to appropriate types
        for key in kwargs:
            try:
                if kwargs[key] == 'True':
                    kwargs[key] = True
                elif kwargs[key] == 'False':
                    kwargs[key] = False
                elif kwargs[key][0].isdigit():
                    if '.' in kwargs[key]:
                        kwargs[key] = float(kwargs[key])
                    elif ',' in kwargs[key] and "label" not in key:
                        kwargs[key] = [int(x) for x in kwargs[key].split(',')]
                    else:
                        kwargs[key] = int(kwargs[key])
                elif kwargs[key] == 'None':
                    kwargs[key] = None
                else:
                    received_data = urllib.parse.unquote(kwargs[key])
                    print("received_data", received_data)
                    kwargs[key] = json.loads(received_data)
            except Exception as e1:
                try:
                    print("key", kwargs[key])
                    received_data = urllib.parse.unquote(kwargs[key])
                    print("received_data", received_data)
                    kwargs[key] = json.loads(received_data)
                except Exception as e2:
                    print("this error is", e1, e2)
                    pass
        
        print(kwargs)

        kwargs['output_type'] = file_type

        print("going to function")
        result = vis_functions[function_name](**kwargs)
        if file_type == 'png':
            
            # result is numpy array png, convert to base64
            if result.max() <= 1:
                result = (result*255).astype(np.uint8)
            img = Image.fromarray(result)
            buffer = BytesIO()
            img.save(buffer, format='PNG')
            buffer.seek(0)
            return send_file(buffer, mimetype='image/png')
        elif file_type == 'svg':
            return result
        else:
            return jsonify({'error': 'File type not supported'}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 400
    

if __name__ == '__main__':
    dash_app.run(debug=True)