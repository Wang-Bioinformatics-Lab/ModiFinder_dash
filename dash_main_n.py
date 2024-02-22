from dash import Dash, html, dcc, Input, Output, State
from Dash_interface import (
    chart_section_n,
    input_section_n,
    url_manager_n,
    computation_n,
)
import sys
from SmallMol_Mod_Site_Localization import ModificationSiteLocator as modSite
from SmallMol_Mod_Site_Localization import Compound_n as compound
from SmallMol_Mod_Site_Localization import visualizer as vis
from SmallMol_Mod_Site_Localization import utils_n as utils
from SmallMol_Mod_Site_Localization import handle_network as hn
import dash_bootstrap_components as dbc

# app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])


if __name__ == "__main__":
    layout = [
        dcc.Store(id="siteLocatorObj"),
        dcc.Store(id="peaksObj"),
        dcc.Store(id="fragmentsObj"),
        dcc.Store(id="InputData"),
        dcc.Store(id="urlData"),
    ]
    args = {
        'USI1': "",
        'USI2': "",
        'SMILES1': "", 
        'SMILES2': "",
        'ppm': 40,
    }

    layout.append(url_manager_n.get_layout())
    layout.append(input_section_n.get_layout(args))
    layout.append(chart_section_n.get_layout())
    app.layout = html.Div(layout)

    url_manager_n.get_callbacks(app, args)
    input_section_n.get_callbacks(app)
    computation_n.get_callbacks(app, modSite, compound, hn)
    chart_section_n.get_callbacks(app, vis, utils)

    app.run_server(debug=True, port=5000, host="0.0.0.0")

server = app.server
