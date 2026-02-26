import traceback

from dash import Dash, html, dcc, Input, Output, State, dash_table
import base64
import pickle
import copy
from urllib.parse import quote

from modifinder import ModiFinder, Compound

from app_utils import get_data, load_helpers, filter_peaks_by_ratio_to_base_peak, adduct_mapping

def get_callbacks(app):
    
    @app.callback(
        [
            Output('siteLocatorObj', 'data'),
            Output('siriusData', 'children'),
            Output('peaksObj', 'data'),
            Output('fragmentsObj', 'data')
         ],
         Output('error-input', 'children'),
        Input('InputData', 'data'),
        )
    def calculate_module(data):
        if data is None:
            return None, None, None, None, None
        
        if data['USI1'] == "" or data['USI2'] == "":
            return None, None, None, None, None
        
        args = copy.deepcopy(data)
        # remove SMILES and USI from args
        args.pop('SMILES1', None)
        args.pop('SMILES2', None)
        usi1 = args.pop('USI1', None)
        usi2 = args.pop('USI2', None)
        
        spectrum1 = get_data(usi1)
        spectrum2 = get_data(usi2)
        if spectrum1.get('adduct') is None:
            # Replace with adduct from data
            spectrum1['adduct'] = adduct_mapping[data['adduct']]    # Should raise error here if we don't know what it is
        if spectrum2.get('adduct') is None:
            # Replace with adduct from data
            spectrum2['adduct'] = adduct_mapping[data['adduct']]

        # TODO: What to do if adduct differs at this point?

        # TODO: Filter adducts in Helpers?

        # Options propagated out of ModiFinder 
        ratio_to_base_peak = args.pop('filter_peaks_variable', None)
        
        # Args to pass to ModiFinder
        args['ppm_tolerance'] = float(args['ppm_tolerance'])
        helper_compounds = args.pop('Helpers', "").strip(' \t\n\r')
        helper_compounds = helper_compounds.replace(" ", "")
        helper_compounds = helper_compounds.split(",")
        # remove empty strings
        helper_compounds = list(filter(None, helper_compounds))
        # remove "" strings
        helper_compounds = list(filter(lambda x: x != "", helper_compounds))
        helper_compounds = load_helpers(
            helper_compounds,
            ratio_to_base_peak=ratio_to_base_peak,
        )

        if data["SMILES1"] == "" or data["SMILES1"] is None:
            data["SMILES1"] = None
        
        if data["SMILES2"] == "" or data["SMILES2"] is None:
            data["SMILES2"] = None

        try:
            # Use known compound adduct
            args['adduct'] = spectrum1.get('adduct', None)

            spectrum1_peaks = spectrum1['peaks']
            spectrum2_peaks = spectrum2['peaks']

            if ratio_to_base_peak:
                spectrum1_peaks = filter_peaks_by_ratio_to_base_peak(spectrum1_peaks, ratio_to_base_peak=ratio_to_base_peak)
                spectrum2_peaks = filter_peaks_by_ratio_to_base_peak(spectrum2_peaks, ratio_to_base_peak=ratio_to_base_peak)

            main_compound = Compound(
                spectrum=spectrum1_peaks,
                precursor_mz=spectrum1['precursor_mz'],
                precursor_charge=spectrum1['precursor_charge'],
                adduct=spectrum1['adduct'],
                smiles=data["SMILES1"]
            )
            mod_compound = Compound(
                spectrum=spectrum2_peaks,
                precursor_mz=spectrum2['precursor_mz'],
                precursor_charge=spectrum2['precursor_charge'],
                adduct=spectrum2['adduct'],
                smiles=data["SMILES2"] if data["SMILES2"] is not None and data["SMILES2"] != "" else None
            )
            
        except Exception as e:
            # if exception is of type value error, return the error message
            if type(e) == ValueError:
                return None, None, None, None, str(e)
            # if exception is of type other, return the error message
            else:
                return None, None, None, None, "Error loading compounds"
        
        if main_compound.structure is None:
            return None, None, None, None, "Error loading SMILES1"

        siteLocator = ModiFinder(main_compound, mod_compound, helpers=helper_compounds, **args)

        peaksObj, fragmentsObj = siteLocator.get_result()
        
        siriusText = "SIRIUS data was not available"
       
        args = copy.deepcopy(data)
        # remove SMILES and USI from args
        args.pop('SMILES1', None)
        args.pop('SMILES2', None)
        args.pop('USI1', None)
        args.pop('USI2', None)

        peaksObj.update({"args": args})

        return base64.b64encode(pickle.dumps(siteLocator)).decode(),  siriusText, base64.b64encode(pickle.dumps(peaksObj)).decode(), base64.b64encode(pickle.dumps(fragmentsObj)).decode(), None