from dash import Dash, html, dcc, Input, Output, State, dash_table, Patch
from dash.exceptions import PreventUpdate
import pickle
import base64
import plotly.graph_objects as go
import rdkit.Chem as Chem
import dash_bootstrap_components as dbc
import numpy as np
from Dash_interface.fragment_selection import FragmentsDisplayAIO
import modifinder.utilities.visualizer as vis
import modifinder.utilities.mol_utils as utils

def dash_svg(text):
    """
    Generates a svg link for dash from a avg text.
    """
    try:
        svg = base64.b64encode(text.encode("utf-8"))
        return "data:image/svg+xml;base64,{}".format(svg.decode())
    except:
        return None

def get_layout():
    # print("in get layout")
    return html.Div(
        id="results",
        children=[
            html.Div(
                children=[
                    dcc.Loading(
                        id="loading-stats",
                        type="default",
                        children=[
                            html.Div(
                                children=[
                                    html.Div(id="stats"),
                                    html.P(id="siriusData"),
                                ],
                                style={
                                    "display": "flex",
                                    "flexDirection": "column",
                                    "justifyContent": "center",
                                    "alignItems": "center",
                                },
                            )
                        ],
                    ),
                    dcc.Loading(
                        id="loading-draw_molecule_result",
                        type="default",
                        children=[
                            html.Div(
                                id="draw_molecule_result",
                                style={
                                    "display": "flex",
                                    "flexDirection": "row",
                                    "justifyContent": "center",
                                    "alignItems": "center",
                                },
                            )
                        ],
                    ),
                ],
                style={
                    "display": "flex",
                    "flexDirection": "row",
                    "justifyContent": "space-around",
                    "alignItems": "center",
                    "position": "relative",
                    "flexWrap": "wrap",
                    "zIndex": "2",
                },
            ),
            dcc.Loading(
                id="loading-peaks",
                type="default",
                children=[
                    html.Div(
                        id="peaks-component",
                        children=[
                            dcc.Graph(id="peaks", style={"width": "100%"}),
                            # html.Div(children = [
                            #     dcc.Slider(
                            #         id='peak_filter_slider',
                            #         min=0,
                            #         max=100,
                            #         step=1,
                            #         value=50,
                            #         marks=None,
                            #         tooltip={"placement": "bottom", "always_visible": True},
                            #     )
                            # ], style = {'width': '30%', 'height': '5vh'}),
                        ],
                        style={"display": "none"},
                    )
                ],
            ),
            dcc.Loading(
                id="loading-peak_info",
                type="default",
                children=[html.Div(id="peak_info")],
            ),
        ],
    )


def get_callbacks(app, diff_to_formula):
    colorsInxMain = {
        "matched_shifted": "red",
        "matched_unshifted": "blue",
        "unmatched": "grey",
        "selected": "green",
    }

    @app.callback(
        [Output("draw_molecule_result", "children"), Output("stats", "children")],
        Input("siteLocatorObj", "data"),
        State("InputData", "data"),
    )
    def update_stats(siteLocatorObj, args):
        if siteLocatorObj == None:
            return None, None
        # print("in update stats")
        siteLocator = pickle.loads(base64.b64decode(siteLocatorObj))

        modified_compound_id = siteLocator._get_unknown()
        
        main_compound_id = siteLocator._get_known_neighbor(modified_compound_id)
        
        main_compound = siteLocator.network.nodes[main_compound_id]['compound']
        modified_compound = siteLocator.network.nodes[modified_compound_id]['compound']

        delta_weight = abs(
                            main_compound.spectrum.precursor_mz
                            - modified_compound.spectrum.precursor_mz
                        )
        ppm_tolerance = siteLocator.ppm_tolerance
        ppm_diff = max(main_compound.spectrum.precursor_mz, modified_compound.spectrum.precursor_mz) * ppm_tolerance / 1000000
        lower_bound = delta_weight - ppm_diff
        upper_bound = delta_weight + ppm_diff
        formulas = diff_to_formula[(diff_to_formula["difference"] >= lower_bound) & (diff_to_formula["difference"] <= upper_bound)]['diff_formula']
        # get unique formulas
        formulas = formulas.drop_duplicates()
        # make comma separated string
        formulas = formulas.str.cat(sep=', ')



        scores = siteLocator.generate_probabilities()
        stats = []
        matched_peaks = siteLocator.get_edge_detail(main_compound.id, modified_compound.id)
        stats.append(("number of matched peaks", str(len(matched_peaks.get_matches_pairs()))))
        # stats.append(("number of shifted", str(len(siteLocator.shifted))))
        stats.append(
            (
                "Δ weight",
                str(round(delta_weight, 5)),
            )
        )
        isAddition = True
        if (
            main_compound.spectrum.precursor_mz
            > modified_compound.spectrum.precursor_mz
        ):
            isAddition = False
        stats.append(("Modification Type", "Addition" if isAddition else "Removal"))
        stats.append(("Possible Modification Formulas", formulas))

        # print("done with stats")
        draw_output = []
        if modified_compound.structure is not None:
            svg1 = vis.draw_modifications(
                modified_compound.structure,
                main_compound.structure,
                output_type="svg",
                highlight_common=False,
                size=(1200, 1200),
                legend_position=(850, 1050),
                legend_font=40,
                modification_categories=['Added', 'Removed'],
            )
            if isAddition:
                trueSite = utils._calculateModificationSites(
                    modified_compound.structure,
                    main_compound.structure,
                    False,
                )[0]
            else:
                trueSite = utils._calculateModificationSites(
                    main_compound.structure,
                    modified_compound.structure,
                    False,
                )[0]
            # stats.append(
            #     (
            #         "Evaluation Score",
            #         str(
            #             round(
            #                 siteLocator.calculate_score(
            #                     trueSite, "average_dist_normalized", scores
            #                 ),
            #                 4,
            #             )
            #         ),
            #     )
            # )
            draw_output.append(
                html.Div(
                    id="prediction_section",
                    children=[
                        html.H3("Modification Highlight "),
                        html.Img(src=dash_svg(svg1), style={"max-width": "30vw"}),
                    ],
                    style={
                        "margin-right": "1vw",
                        "padding-right": "1vw",
                        "border-right": "1px dashed #aaa",
                    },
                )
            )


        svg2 = vis.draw_molecule_heatmap(main_compound.structure, scores, output_type="svg", shrink_labels=True, show_labels=True, size=(1200, 1200))
        # print("making svg3")
        # svg3 = ""
        svg3 = vis.draw_molecule_heatmap(main_compound.structure, scores, output_type="svg", shrink_labels=True, show_labels=True, show_legend=False, size=(1200, 1200))
        # print("done with svg3")
        
        draw_output.append(
            html.Div(
                id="substr_section",
                children=[
                    html.H3("Prediction"),
                    html.Img(src=dash_svg(svg2), style={"max-width": "30vw"}),
                    # download button to download the svg
                    html.Div(
                        [
                            dbc.Button(
                                [html.I(className='bi bi-download', style={'font-size': '1rem'})], # House icon
                                title="Download SVG",
                                id="download_svg_button",
                                color="secondary",
                                style={"position": "absolute", "top": "5%", "right": "5%"},
                            ),
                            # store the svg in a hidden div
                            html.Div(id="hidden_svg", children=svg3, style={"display": "none"}),
                            dcc.Download(id="download_heatmap_svg"),
                        ]
                    ),
                ],
                style={
                    "position": "relative",
                }
            )
        )

        stats_section = []
        for stat in stats:
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

        stats = dbc.Card(
            children=stats_section,
            style={"flex": "1", "width": "100%", "padding": "10px"},
        )

        return draw_output, stats

    @app.callback(
        [Output("peaks", "figure"), Output("peaks-component", "style")],
        Input("peaksObj", "data"),
        # Input('peak_filter_slider', 'value'),
        prevent_initial_call=True,
    )
    def update_peaks(data):  # , slider_value):
        if data == None:
            return {}, {"display": "none"}
        peaksObj = pickle.loads(base64.b64decode(data))
        main_compound_peaks = peaksObj["main_compound_peaks"]
        mod_compound_peaks = peaksObj["mod_compound_peaks"]
        matched_peaks = peaksObj["matched_peaks"]
        args = peaksObj["args"]
        main_precursor_mz = peaksObj["main_precursor_mz"]
        mod_precursor_mz = peaksObj["mod_precursor_mz"]

        fig = go.Figure()
        typesInxMain = {"matched_shifted": [], "matched_unshifted": [], "unmatched": []}
        x1 = []
        y1 = []
        for peak in main_compound_peaks:
            x1.append(peak[0])
            y1.append(peak[1])

        # topPeakCount = slider_value
        topPeakCount = max(
            len(main_compound_peaks),
            len(mod_compound_peaks),
        )
        topPeaksInxModif = sorted(range(len(y1)), key=lambda i: y1[i])[-topPeakCount:]
        hoverData = {"main": [], "modified": []}
        for i in topPeaksInxModif:
            flag = False
            for j in matched_peaks:
                if j[0] == i:
                    if (
                        abs(
                            main_compound_peaks[i][0]
                            - mod_compound_peaks[j[1]][0]
                        )
                        > args["mz_tolerance"]
                    ):
                        typesInxMain["matched_shifted"].append(i)
                        hoverData["main"].append(j[1])
                    else:
                        typesInxMain["matched_unshifted"].append(i)
                    flag = True
                    break
            if not flag:
                typesInxMain["unmatched"].append(i)

        typesInxModified = {
            "matched_shifted": [],
            "matched_unshifted": [],
            "unmatched": [],
        }
        x2 = []
        y2 = []
        for peak in mod_compound_peaks:
            x2.append(peak[0])
            y2.append(peak[1])

        topPeaksInxModif = sorted(range(len(y2)), key=lambda i: y2[i])[-topPeakCount:]
        for i in topPeaksInxModif:
            flag = False
            for j in matched_peaks:
                if j[1] == i:
                    if (
                        abs(
                            main_compound_peaks[j[0]][0]
                            - mod_compound_peaks[j[1]][0]
                        )
                        > 0.1
                    ):
                        typesInxModified["matched_shifted"].append([i, j[0]])
                        hoverData["modified"].append(j[0])
                    else:
                        typesInxModified["matched_unshifted"].append([i, j[0]])
                    flag = True
                    break
            if not flag:
                typesInxModified["unmatched"].append([i, -1])

        minX = min(min(x1), min(x2))
        maxX = max(max(x1), max(x2))
        minX = min(minX, main_precursor_mz, mod_precursor_mz)
        maxX = max(maxX, main_precursor_mz, mod_precursor_mz)

        for inx_type in typesInxMain:
            x_main = [round(x1[j], 4) for j in typesInxMain[inx_type]]
            y1_ = [y1[j] for j in typesInxMain[inx_type]]
            y_main = [y / max(y1_) * 100 for y in y1_]
            x_modified = [round(x2[j[0]], 4) for j in typesInxModified[inx_type]]
            y2_ = [y2[j[0]] for j in typesInxModified[inx_type]]
            y_modified = [-j / max(y2_) * 100 for j in y2_]
            indicis = typesInxMain[inx_type] + [
                j[0] for j in typesInxModified[inx_type]
            ]
            x_ = x_main + x_modified
            y_ = y_main + y_modified
            colors = [colorsInxMain[inx_type]] * len(x_)
            if inx_type == "unmatched":
                visibility = "legendonly"
                if len(typesInxModified["matched_shifted"]) == 0 and len(
                    typesInxModified["matched_unshifted"]
                ) == 0:
                    visibility = None
                
                fig.add_trace(
                    go.Bar(
                        x=x_,
                        y=y_,
                        width=(maxX - minX) / 500,
                        hovertext=indicis,
                        name=inx_type,
                        visible=visibility,
                        marker_color=colors,
                    )
                )
            elif inx_type == "matched_shifted":
                hovertext = []
                for i in range(len(x_main)):
                    hovertext.append(
                        str(indicis[i])
                        + " "
                        + "matched to:"
                        + str(hoverData["main"][i])
                    )
                for i in range(len(x_main), len(x_main) + len(x_modified)):
                    hovertext.append(
                        str(indicis[i])
                        + " "
                        + "matched to:"
                        + str(hoverData["modified"][i - len(x_main)])
                    )

                fig.add_trace(
                    go.Bar(
                        x=x_,
                        y=y_,
                        hovertext=hovertext,
                        name=inx_type,
                        width=(maxX - minX) / 500,
                        marker_color=colors,
                    )
                )
            else:
                fig.add_trace(
                    go.Bar(
                        x=x_,
                        y=y_,
                        hovertext=indicis,
                        name=inx_type,
                        width=(maxX - minX) / 500,
                        marker_color=colors,
                    )
                )

        # add vertical dashed line for the precursor m/z
        fig.add_trace(
            go.Scatter(
                x=[main_precursor_mz, main_precursor_mz],
                y=[0, 100],
                mode="lines",
                line=go.scatter.Line(color="black", dash="dash", width= (maxX - minX) / 600),
                name='known precursor m/z',
                # showlegend=False,
                # hoverinfo='skip'
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[mod_precursor_mz, mod_precursor_mz],
                y=[-100, 0],
                mode="lines",
                line=go.scatter.Line(color="black", dash="dot", width= (maxX - minX) / 600),
                name='modified precursor m/z',
                # showlegend=False,
                # hoverinfo='skip'
            )
        )

        # minX = min(minX, main_precursor_mz, mod_precursor_mz)
        # maxX = max(maxX, main_precursor_mz, mod_precursor_mz)

        # fig.update_traces(
        #     width=(maxX - minX) / 400,
        # )
        fig.update_layout(
            title="Alignment of Peaks",
            bargap=0,
            xaxis_title="m/z",
            yaxis_title="Intensity",
            legend_title="Peak Type",
            # remove paddings
            margin=dict(l=0, r=0, t=40, b=0),

            # range of xaxis
            xaxis=dict(
                range=[minX - (maxX * 0.01), maxX * 1.01],
            ),
        )
        return fig, {
            "display": "flex",
            "position": "relative",
            "flexDirection": "column",
            "justifyContent": "center",
            "width": "98%",
            "alignItems": "center",
            "padding": "20px",
            "margin": "auto",
            "margin-top": "1vh",
            "zIndex": "1",
        }


    # @app.callback(
    #     Output("peak_info", "children", allow_duplicate=True),
    #     Input("siteLocatorObj", "data"), 
    #     prevent_initial_call=True,
    # )
    # def clear_peak_info(data):
    #     if data == None:
    #         return ""
    #     else:
    #         return "Select a peak to see its fragments"

    @app.callback(
        Output("peak_info", "children", allow_duplicate=True),
        Input("peaks", "clickData"),
        State("fragmentsObj", "data"),
        prevent_initial_call=True,
    )
    def display_click_data(clickData, fragmentsObj):
        if clickData:
            try:
                fragmentsObj = pickle.loads(base64.b64decode(fragmentsObj))

                clicked_peak_x = float(clickData["points"][0]["x"])
                clicked_peak_y = float(clickData["points"][0]["y"])
                if clicked_peak_y < 0:
                    return "only known compound peaks are clickable"
                if clicked_peak_x > fragmentsObj["Precursor_MZ"] - 1:
                    return "precursor peak, " + str(clicked_peak_x)
                
                structure = fragmentsObj["structure"]
                frags_map = fragmentsObj["frags_map"]
                peaks = fragmentsObj["peaks"]

                peak_index = -1
                for i, peak in enumerate(peaks):
                    if abs(peak[0]- clicked_peak_x)/clicked_peak_x*1000000 < 40:
                        peak_index = i
                        break
                if peak_index == -1:
                    return "error in finding peak index"
                
                fragments = list(frags_map[peak_index])
                result_posibility_indicies = []
                for fragment in fragments:
                    fragment_indicies = []
                    for i in range(len(structure.GetAtoms())):
                        if fragment & 1 << i:
                            fragment_indicies.append(i)
                    result_posibility_indicies.append(fragment_indicies)
                
                return FragmentsDisplayAIO(
                    result_posibility_indicies,
                    structure,
                    {
                        "mz": clickData["points"][0]["x"],
                        "intensity": clickData["points"][0]["y"],
                        "all_fragments": fragments,
                    },
                    "fragmentDisplay",
                )
            except:
                import traceback

                traceback.print_exc()
                return "siteLocator object not found"
        return None

    # change the color of the bar when clicked
    @app.callback(
        Output("peaks", "figure", allow_duplicate=True),
        Input("peaks", "clickData"),
        State("peaks", "figure"),
        prevent_initial_call=True,
    )
    def change_bar_color(clickData, figure):
        # set all colors to default
        patched_figure = Patch()
        for i in range(len(figure["data"])):
            for j in range(len(figure["data"][i]["x"])):
                if figure["data"][i]["name"] in colorsInxMain and figure["data"][i]["marker"]["color"][j] != colorsInxMain[figure["data"][i]["name"]]:
                    patched_figure["data"][i]["marker"]["color"][j] = colorsInxMain[figure["data"][i]["name"]]
        
        if clickData:
            try:
                clicked_peak_x = float(clickData["points"][0]["x"])
                clicked_peak_y = float(clickData["points"][0]["y"])
                if clicked_peak_y < 0:
                    return figure
                for i in range(len(figure["data"])):
                    for j in range(len(figure["data"][i]["x"])):
                        if (figure["data"][i]["x"][j] == clicked_peak_x) and (
                            figure["data"][i]["y"][j] == clicked_peak_y
                        ):
                            patched_figure["data"][i]["marker"]["color"][j] = "green"
                            # figure["data"][i]["marker"]["color"][j] = "green"
                            # if matched shifted peak, highlight the corresponding peak in the other bar
                            if figure["data"][i]["name"] == "matched_shifted":
                                index = figure["data"][i]["hovertext"][j].split(":")[1]
                                for l in range(len(figure["data"][i]["x"])):
                                    if (figure["data"][i]["hovertext"][l].split(" ")[0] == index and figure["data"][i]["y"][l] < 0):
                                        patched_figure["data"][i]["marker"]["color"][l] = "olive"
                                        break

                            break
            except:
                import traceback

                traceback.print_exc()
                return patched_figure

        return patched_figure

    @app.callback(
        [Output("siteLocatorObj", "data", allow_duplicate=True), 
         Output("peak_info", "children", allow_duplicate=True),
         Output('fragmentsObj', 'data', allow_duplicate=True)],
        Input(FragmentsDisplayAIO.ids.fragment_data("fragmentDisplay"), "data"),
        State("siteLocatorObj", "data"),
        prevent_initial_call=True,
    )
    def apply_structure_filter(data, siteLocatorObj):
        if siteLocatorObj == None:
            raise PreventUpdate
        if data["has_changed"] == False:
            raise PreventUpdate
        siteLocator = pickle.loads(base64.b64decode(siteLocatorObj))
        
        modified_compound_id = siteLocator._get_unknown()
        main_compound_id = siteLocator._get_known_neighbor(modified_compound_id)
        main_compound = siteLocator.network.nodes[main_compound_id]['compound']
        main_compound_peaks = [(main_compound.spectrum.mz[i], main_compound.spectrum.intensity[i]) for i in range(len(main_compound.spectrum.mz))]
        modified_compound = siteLocator.network.nodes[modified_compound_id]['compound']
        
        ind = main_compound.spectrum.get_peak_indexes(data["mz"])
        main_compound.peak_fragments_map[ind[0]] = [data["all_fragments"][i] for i in data["selected_fragments"]]
        
        fragmentsObj = {
            "frags_map": main_compound.peak_fragments_map,
            "structure": main_compound.structure,
            "peaks": main_compound_peaks,
            "Precursor_MZ": main_compound.spectrum.precursor_mz,
        }


        fragments = list(main_compound.peak_fragments_map[ind[0]])
        result_posibility_indicies = []
        for fragment in fragments:
            fragment_indicies = []
            for i in range(len(main_compound.structure.GetAtoms())):
                if fragment & 1 << i:
                    fragment_indicies.append(i)
            result_posibility_indicies.append(fragment_indicies)

        new_peak_info = FragmentsDisplayAIO(
            result_posibility_indicies,
            main_compound.structure,
            {
                "mz": data["mz"],
                "intensity": data["intensity"],
                "all_fragments": fragments,
            },
            "fragmentDisplay",
        )

        return base64.b64encode(pickle.dumps(siteLocator)).decode(), new_peak_info, base64.b64encode(pickle.dumps(fragmentsObj)).decode()
    
    @app.callback(
        Output("download_heatmap_svg", "data"),
        Input("download_svg_button", "n_clicks"),
        State("hidden_svg", "children"),
    )
    def download_svg(n_clicks, svg):
        if n_clicks:
            return dict(content=svg, filename="heatmap.svg")
        return None
