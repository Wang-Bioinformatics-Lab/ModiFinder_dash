import dash
from dash import Dash, dcc, html, Input, Output, State, MATCH, ALL
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
from modifinder.utilities import spectrum_utils


from flask import Flask, send_file, request, jsonify
import json
from app import app
from furl import furl
from myopic_mces import MCES
supported_adducts =  ['[M+H]1+', '[M-H]1-', '[M+Na]1+', '[M+NH4]1+', '[M+K]1+', '[M+Cl]1-', '[M+Br]1-']


def png_to_showable_src(png):
    if png.max() <= 1:
        png = (png*255).astype(np.uint8)
    img = Image.fromarray(png)
    buffer = BytesIO()
    img.save(buffer, format='PNG')
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    return f'data:image/png;base64,{img_base64}'


def get_spectrum(input_text, adduct):
    try:
        spectrum = Spectrum(input_text)
    except Exception as e:
        try:
            peaks = input_text.split('\n')
            # remove empty lines
            peaks = [peak for peak in peaks if peak]
            
            # remove starting and trailing whitespaces
            peaks = [peak.strip() for peak in peaks]
            
            # replace tabs with spaces
            peaks = [peak.replace('\t', ' ') for peak in peaks if peak]
            
            # remove multiple spaces
            peaks = [peak.split(' ') for peak in peaks]
            for peak in peaks:
                while '' in peak:
                    peak.remove('')
            peaks = [' '.join(peak) for peak in peaks]
            mz = [float(peak.split(' ')[0]) for peak in peaks]
            intensity = [float(peak.split(' ')[1]) for peak in peaks]
            
            spectrum = Spectrum(mz = mz, intensity = intensity, adduct = adduct)
            
        except Exception as e2:
            raise Exception(f"Error: {e} | {e2}")
    
    return spectrum
    

# dash.register_page(__name__)
# app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP])
dash_app = Dash(
    name="Merge dashboard",
    server=app,
    url_base_pathname="/merge/",
    external_stylesheets=[dbc.themes.BOOTSTRAP],
)

dash_app.layout = html.Div([
    html.H1("Merge Spectra Peaks"),
    dbc.Button("Add Spectrum", id="add-spectrum", n_clicks=0, style={'margin': '10px'}),
    html.Div(id="spectra-container", children=[]),
    # add ppm input
    html.Div([
        dbc.InputGroup(
            [dbc.InputGroupText("PPM Tolerance"),
            dbc.Input(placeholder="ppm", type="number", id="ppm-tolerance", value=10)],
            style = {'margin': '10px'}
        ),
        dbc.InputGroup(
            # concensus ratio (0-1)
            [dbc.InputGroupText("Concensus Ratio"),
            dbc.Input(placeholder="ratio", type="number", id="concensus-ratio", value=0.5)],
            style = {'margin': '10px'}
        )
    ], style = {'display': 'flex', 'flexDirection': 'row', 'alignItems': 'center', 'justifyContent': 'center'}),
    dbc.Button("Done", id="done", n_clicks=0, style={'margin': '10px'}),
    
    html.Div([
        dcc.Textarea(id="output", style={"minHeight": "150px", 'minWidth': "250px"}, readOnly=True),
        dcc.Input(id="output_adduct", style={"height": "150px", "border": "none"}, readOnly=True),
        html.Div(id="output_image"),
    ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'center', 'justifyContent': 'center'}),
    
    html.Div([
        html.P("you can use "),
        html.A(" TinyMass ", href="https://www.chemcalc.org/chemcalc", target="_blank"),
        html.P(" and copy the merged peaks calculated above to get a USI"),
    ], style={'width': '80%', 'margin': 'auto', 'textAlign': 'center', 'display': 'flex', 
              'flexDirection': 'row', 'alignItems': 'baseline', 'justifyContent': 'center'}),
], style={'width': '80%', 'margin': 'auto', 'textAlign': 'center'})

@dash_app.callback(
    Output("spectra-container", "children"),
    Input("add-spectrum", "n_clicks"),
    State("spectra-container", "children")
)
def add_spectrum(n_clicks, current_spectra):
    if current_spectra is None:
        current_spectra = []
    # Each spectrum will be identified by its index (n_clicks used as id)
    new_spectrum_id = n_clicks if n_clicks > 0 else 0
    new_spectrum = html.Div([
        dbc.InputGroup(
            [dbc.InputGroupText(f"Spectrum {new_spectrum_id + 1}"),
            #  multiple line input for peaks
             dbc.Textarea(placeholder="peaks", id={"type": "spectrum-name", "spectrum": new_spectrum_id}, style={'minHeight': '100px'})]),
        dbc.InputGroup(
            [dbc.InputGroupText("Adduct"),
             dcc.Dropdown(
               supported_adducts, supported_adducts[0],
                id = {"type": "adduct", "spectrum": new_spectrum_id},
                clearable=False,
                style = {'minWidth': '120px', 'borderTopLeftRadius': 0, 'borderBottomLeftRadius': 0})],
            style = {'margin': '10px'}),
        # add image of spectrum
        html.Div(id={"type": "spectrum-illustration", "spectrum": new_spectrum_id}),
        ],style = {'display': 'flex', 'flexDirection': 'row', 'alignItems': 'center', 'justifyContent': 'center'}
        )
    # new_spectrum = html.Div([
    #     html.H3(f"Spectrum {new_spectrum_id + 1}"),
    #     dcc.Input(type="text", placeholder="peaks", id={"type": "spectrum-name", "spectrum": new_spectrum_id}, style={'margin': '10px'}),
    # ], style={'border': '1px solid #ccc', 'padding': '10px', 'margin': '10px'})
    
    current_spectra.append(new_spectrum)
    return current_spectra

# callback to update the image of the spectrum, if can't parse peaks, show error
@dash_app.callback(
    Output({"type": "spectrum-illustration", "spectrum": MATCH}, "children"),
    # Output({"type": "adduct", "spectrum": MATCH}, "value"),
    Input({"type": "spectrum-name", "spectrum": MATCH}, "value"),
    State({"type": "adduct", "spectrum": MATCH}, "value"),
)
def update_spectrum_image(peaks, adduct):
    if not peaks:
        return ""
    
    spectrum = None
    try:
        spectrum = get_spectrum(peaks, adduct)
    except Exception as e:
        error_item = html.Div([
            html.P(f"{e}", style={'color': 'red'}),
        ])
        return error_item
    
    img = mf_vis.draw_spectrum(spectrum, bar_width = 1)
    img_src = png_to_showable_src(img)
    img_item = html.Img(src=img_src, style={'height': '200px'})
    return img_item
    


# Callback to collect all peak values when "Done" is clicked
@dash_app.callback(
    Output("output", "value"),
    Output("output_adduct", "value"),
    Output("output_image", "children"),
    Input("done", "n_clicks"),
    State({"type": "spectrum-name", "spectrum": ALL}, "value"),
    State({"type": "adduct", "spectrum": ALL}, "value"),
    State("ppm-tolerance", "value"),
    State("concensus-ratio", "value")
)
def show_peaks(n_clicks, spectra, adducts, ppm, concensus_ratio):
    if not n_clicks:
        return ""
    
    spectrums = []
    for i, peaks in enumerate(spectra):
        try:
            spectrum = get_spectrum(peaks, adducts[i])
            spectrums.append(spectrum)
        except Exception as e:
            return f"Error: {e}", None, None

    try:
        merged_spectrum = spectrum_utils.aggregate_spectrums(spectrums, ppm_tolerance = ppm,
                                                            consensus_majority_ratio = concensus_ratio)
    except Exception as e:
        return f"Error: {e}", None, None
    
    # output peaks
    output = []
    for i, mz in enumerate(merged_spectrum.mz):
        output.append(f"{round(mz, 4)} {round(merged_spectrum.intensity[i], 4)}")
    output = '\n'.join(output)
    
    # output adduct
    output_adduct = f"Adduct: {merged_spectrum.adduct}"
    
    # output image
    img = mf_vis.draw_spectrum(merged_spectrum, bar_width = 1)
    img_src = png_to_showable_src(img)
    img_item = html.Img(src=img_src, style={'height': '200px'})
    
    return output, output_adduct, img_item
    
    
    

if __name__ == '__main__':
    dash_app.run_server(debug=True)