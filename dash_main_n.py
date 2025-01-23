import argparse

from dash import Dash, html, dcc, Input, Output, State
from flask import request

from Dash_interface import (
    chart_section_n,
    input_section_n,
    url_manager_n,
    computation_n,
)
import sys
from modifinder import ModiFinder, Compound
import dash_bootstrap_components as dbc
import pandas as pd
import json

from app import app

# app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

#app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP])
dash_app = Dash(
    name="modifinder dashboard",
    server=app,
    url_base_pathname="/",
    external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
)

dash_app.title = "ModiFinder Dashboard"

# setting tracking token
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

ACKNOWLEDGEMENTS = [
    dbc.CardHeader(html.H5("Acknowledgements")),
    dbc.CardBody(
        [
            "This work was supported by Agilent Technologies"
        ]
    )
]

def _setup():
    print("LOADING SETUP in MODIFINDER", file=sys.stderr, flush=True)

    diff_to_formula = pd.read_csv("data/meta/formula_diff.csv")
    helpers = json.load(open("data/meta/helpers.json"))

    layout = [
        NAVBAR,
        dcc.Store(id="siteLocatorObj"),
        dcc.Store(id="peaksObj"),
        dcc.Store(id="fragmentsObj"),
        dcc.Store(id="InputData"),
        dcc.Store(id="urlData"),
    ]
    initial_inputs = {
        'USI1': "",
        'USI2': "",
        'SMILES1': "", 
        'SMILES2': "",
        'Helpers': "",
        'adduct': "[M+H]1+",
        'ppm_tolerance': 40,
        'filter_peaks_variable': 0.01,
    }

    layout.append(url_manager_n.get_layout())
    layout.append(input_section_n.get_layout(initial_inputs))
    layout.append(chart_section_n.get_layout())
    layout.append(dbc.Card(CONTRIBUTORS_DASHBOARD, style={"width": "98%", "margin": "1vh auto"}))
    layout.append(dbc.Card(ACKNOWLEDGEMENTS, style={"width": "98%", "margin": "1vh auto"}))
    dash_app.layout = html.Div(layout)

    url_manager_n.get_callbacks(dash_app, initial_inputs)
    input_section_n.get_callbacks(dash_app, helpers)
    computation_n.get_callbacks(dash_app)
    chart_section_n.get_callbacks(dash_app, diff_to_formula)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Dash server")
    parser.add_argument("--port", type=int, default=5002)
    # add boolean argument to enable/disable debug mode
    parser.add_argument("--debug", action="store_true", default=False)
    args = parser.parse_args()

    print("RUNNING MODIFINDER MAIN", file=sys.stderr, flush=True)

    _setup()

    dash_app.run_server(debug=args.debug, port=args.port, host="0.0.0.0")
else:
    print("RUNNING MODIFINDER SETUP NOT MAIN", file=sys.stderr, flush=True)
    _setup()
