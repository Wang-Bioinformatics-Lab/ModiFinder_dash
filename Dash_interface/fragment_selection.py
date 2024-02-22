from dash import Dash, Output, Input, State, html, dcc, callback, MATCH, no_update
from dash.exceptions import PreventUpdate
import uuid
import base64
from dash.development.base_component import Component
from SmallMol_Mod_Site_Localization import visualizer
import dash_bootstrap_components as dbc


def dash_svg(text):
    """
    Generates a svg link for dash from a avg text.
    """
    try:
        svg = base64.b64encode(text.encode("utf-8"))
        return "data:image/svg+xml;base64,{}".format(svg.decode())
    except:
        return None


# class SingleFragmentAIO(html.Div):  # html.Div will be the "parent" component
#     def __init__(self, fragment, mol, formula = "None", mode="display", aio_id=None,*args, **kwargs):
#         if aio_id is None:
#             aio_id = str(uuid.uuid4())

#         self.image = visualizer.highlightMolIndices(mol, fragment)
#         self.formula = formula

#         super().__init__(children=[
#             html.Img(src= dash_svg(self.image), style={"width": "300px"}),
#             html.Div(self.formula) if self.formula != "None" else None,
#             # show the selected state of the fragment if it is in edit mode
#             dcc.Checklist(
#                 options=[
#                     {'label': 'Selected', 'value': True}
#                 ],
#                 value=[True],
#                 id={
#                     "component": "FragmentsDisplayAIO",
#                     "subcomponent": "SingleFragmentAIO",
#                     "aio_id": aio_id,
#                 }
#             ) if mode == "edit" else None
#         ], *args, **kwargs)

# @property
# def image(self):
#     return visualizer.draw_fragment(self.fragment, self.mol)

# @property
# def formula(self):
#     return self.fragment.formula


# All-in-One Components should be suffixed with 'AIO'
class FragmentsDisplayAIO(html.Div):  # html.Div will be the "parent" component
    class ids:
        button = lambda aio_id: {
            "component": "FragmentsDisplayAIO",
            "subcomponent": "button",
            "aio_id": aio_id,
        }
        fragment_selection = lambda aio_id: {
            "component": "FragmentsDisplayAIO",
            "subcomponent": "fragment-selection",
            "aio_id": aio_id,
        }
        fragment_data = lambda aio_id: {
            "component": "FragmentsDisplayAIO",
            "subcomponent": "store",
            "aio_id": aio_id,
        }
        meta_text = lambda aio_id: {
            "component": "FragmentsDisplayAIO",
            "subcomponent": "meta-text",
            "aio_id": aio_id,
        }
        all_handler = lambda aio_id: {
            "component": "FragmentsDisplayAIO",
            "subcomponent": "all-handler",
            "aio_id": aio_id,
        }
        select_all = lambda aio_id: {
            "component": "FragmentsDisplayAIO",
            "subcomponent": "select-all",
            "aio_id": aio_id,
        }
        deselect_all = lambda aio_id: {
            "component": "FragmentsDisplayAIO",
            "subcomponent": "deselect-all",
            "aio_id": aio_id,
        }

    # Make the ids class a public class
    ids = ids

    def __init__(self, fragments_indicies, mol, info, aio_id=None, *args, **kwargs):
        if aio_id is None:
            aio_id = str(uuid.uuid4())

        # self.mode = "display"
        # self.mol = mol
        # self.fragments = fragments
        # self.selected_fragments = [i for i in range(len(fragments))]

        data = {
            "mz": info["mz"],
            "intensity": info["intensity"],
            "all_fragments": info["all_fragments"],
            "selected_fragments": [i for i in range(len(fragments_indicies))],
            "has_changed": False,
        }

        infoSection = html.Div(
            id="peak_stats",
            children=[
                html.H3("Peak Info"),
                html.Div(
                    children=[
                        html.P("m/z: " + str(info["mz"])),
                        html.P("Intensity: " + str(info["intensity"])),
                    ]
                ),
            ],
        )

        self.aio_id = aio_id
        super().__init__(
            children=[
                # dbc.Button(self.mode, id= self.ids.button(self.aio_id)),
                infoSection,
                # show button to edit the fragment selection if there is fragment
                dbc.Button("click to edit", id=self.ids.button(self.aio_id)) if len(fragments_indicies) > 0 else html.P("No fragments found"),
                html.Div(id= self.ids.all_handler(self.aio_id), children=[
                    dbc.Button("Select all", id= self.ids.select_all(self.aio_id), n_clicks=0),
                    dbc.Button("Deselect all", id= self.ids.deselect_all(self.aio_id), n_clicks=0),
                ], style={"display": "none"}),

                dcc.Checklist(
                    [
                        {
                            "label": html.Img(
                                src=dash_svg(
                                    visualizer.highlightMolIndices(mol, fragment)
                                ),
                                style={"width": "300px"},
                            ),
                            "value": i,
                            "disabled": True,
                        }
                        for i, fragment in enumerate(fragments_indicies)
                    ],
                    value=[i for i in range(len(fragments_indicies))],
                    id=self.ids.fragment_selection(self.aio_id),
                    inline=True,
                ),
                dcc.Store(id=self.ids.fragment_data(self.aio_id), data=data),
            ], style={"margin-bottom": "20px"}
        )

    # deselect all fragments
    @callback(
        Output(ids.fragment_selection(MATCH), "value", allow_duplicate=True),
        Input(ids.deselect_all(MATCH), "n_clicks"),
        prevent_initial_call=True,
    )
    def deselect_all(n_clicks):
        if n_clicks is None or n_clicks < 1:
            raise PreventUpdate
        return []
    
    # select all fragments
    @callback(
        Output(ids.fragment_selection(MATCH), "value"),
        Input(ids.select_all(MATCH), "n_clicks"),
        State(ids.fragment_data(MATCH), "data"),
        prevent_initial_call=True,
    )
    def select_all(n_clicks, data):
        if n_clicks is None or n_clicks < 1:
            raise PreventUpdate
        return [i for i in range(len(data["all_fragments"]))]

    # change the mode of the fragment selection
    @callback(
        Output(ids.button(MATCH), "children"),
        Output(ids.fragment_selection(MATCH), "options"),
        Output(ids.fragment_data(MATCH), "data"),
        Output(ids.all_handler(MATCH), "style"),
        Input(ids.button(MATCH), "n_clicks"),
        State(ids.button(MATCH), "children"),
        State(ids.fragment_selection(MATCH), "options"),
        State(ids.fragment_data(MATCH), "data"),
        State(ids.fragment_selection(MATCH), "value"),
        prevent_initial_call=True,
    )
    def change_mode(n_clicks, mode, options, data, selected_fragments):
        # print("change_mode", n_clicks)
        if n_clicks is None or n_clicks < 1:
            raise PreventUpdate
        if mode == "click to edit":
            mode = "click to apply"
            options = [
                {"label": option["label"], "value": option["value"]}
                for option in options
            ]
        else:
            mode = "click to edit"
            options = [
                {"label": option["label"], "value": option["value"], "disabled": True}
                for option in options
            ]
            # check if the selected fragments have changed
            if data["selected_fragments"] != selected_fragments:
                data["selected_fragments"] = selected_fragments
                data["has_changed"] = True
                # print("selected_fragments changed")
            # else:
                # print("selected_fragments not changed")
            # sort to show the selected fragments first
            options.sort(key=lambda option: option["value"] not in selected_fragments)
        
        if mode == "click to edit":
            style = {"display": "none"}
        else:
            style = {"display": "block"}
        return mode, options, data, style
