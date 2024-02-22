from dash import Dash, html, dcc, Input, Output, State, dash_table
import base64
import pickle
import json

def get_callbacks(app, modSite, Compound, hn):
    
    @app.callback(
        [Output('siteLocatorObj', 'data'),  Output('siriusData', 'children'), Output('peaksObj', 'data'), Output('fragmentsObj', 'data')],
        Input('InputData', 'data'),
        )
    def calculate_module(data):
        if data is None:
            return None, None, None, None
        
        if data['USI1'] == "" or data['USI2'] == "":
            return None, "Please input USI", None, None
        try:
            main_info = hn.getDataFromUsi(data['USI1'])
            main_info['Adduct'] = "M+H"
            mod_info = hn.getDataFromUsi(data['USI2'])
            mod_info['Adduct'] = "M+H"
            # print("data loaded")
        except:
            return None, "error loading data for USI", None, None
        
        
        main_compound = Compound.Compound(main_info, data['SMILES1'], data.get('args', {}))
        sirius_data = None
        try:
            specID = data['USI1'].split(":")[-1]
            with open("./data/SIRIUS/"+specID+"_fragmentationtree.json", "rb") as f:
                sirius_data = json.load(f)
            main_compound.apply_sirius(sirius_data)
        except:
            pass

        main_compound.fragments.fragment_masses = []
        main_compound.fragments.fragment_info = []

        # print("main compound loaded")
        mod_compound_args = data.get('args', {})
        mod_compound_args['should_fragment'] = False
        mod_compound = Compound.Compound(mod_info, data.get('SMILES2', None), args = mod_compound_args)
        # print("mod compound loaded")
        siteLocator = modSite.ModificationSiteLocator(main_compound, mod_compound, data.get('args', {}))
        # print("site locator loaded")
        siriusText = "SIRIUS data was not available"
        if sirius_data is not None:
            siriusText = "SIRIUS data was available"
        # else:
        #     print("SIRIUS data was not available", data['USI1'])
        # if siteLocator.main_compound.Precursor_MZ > siteLocator.modified_compound.Precursor_MZ:
        #     return None, "Molecule precursor mass is higher than modified precursor mass", siriusText
        # else:
        peaksObj = {
            "main_compound_peaks": siteLocator.main_compound.peaks,
            "mod_compound_peaks": siteLocator.modified_compound.peaks,
            "matched_peaks": siteLocator.matched_peaks,
            "args": data.get('args', {}),
            "main_precursor_mz": siteLocator.main_compound.Precursor_MZ,
            "mod_precursor_mz": siteLocator.modified_compound.Precursor_MZ,
        }

        fragmentsObj = {
            "frags_map": siteLocator.main_compound.peak_fragments_map,
            "structure": siteLocator.main_compound.structure,
            "peaks": siteLocator.main_compound.peaks,
            "Precursor_MZ": siteLocator.main_compound.Precursor_MZ,
        }

        return base64.b64encode(pickle.dumps(siteLocator)).decode(),  siriusText, base64.b64encode(pickle.dumps(peaksObj)).decode(), base64.b64encode(pickle.dumps(fragmentsObj)).decode()