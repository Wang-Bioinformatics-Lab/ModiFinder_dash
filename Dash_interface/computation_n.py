from dash import Dash, html, dcc, Input, Output, State, dash_table
import base64
import pickle
import json
import copy
from modifinder import ModiFinder, Compound
from rdkit import Chem

def get_callbacks(app):
    
    @app.callback(
        [Output('siteLocatorObj', 'data'),  Output('siriusData', 'children'), Output('peaksObj', 'data'), Output('fragmentsObj', 'data')], Output('error-input', 'children'),
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
        args.pop('USI1', None)
        args.pop('USI2', None)
        args["normalize_peaks"] = True
        args["remove_large_peaks"] = True
        args["ratio_to_base_peak"] = float(args["filter_peaks_variable"])
        args['ppm_tolerance'] = float(args['ppm_tolerance'])
        helper_compounds = args.pop('Helpers', "").strip(' \t\n\r')
        helper_compounds = helper_compounds.replace(" ", "")
        helper_compounds = helper_compounds.split(",")
        # remove empty strings
        helper_compounds = list(filter(None, helper_compounds))
        # remove "" strings
        helper_compounds = list(filter(lambda x: x != "", helper_compounds))

        if data["SMILES1"] == "" or data["SMILES1"] is None:
            data["SMILES1"] = None
        
        if data["SMILES2"] == "" or data["SMILES2"] is None:
            data["SMILES2"] = None

        try:
            if data['adduct']:
                args['adduct'] = data['adduct']
            main_compound = Compound(data['USI1'], **args)
            if data["SMILES1"] is not None:
                main_compound.update(smiles=data["SMILES1"])
            mod_compound = Compound(data['USI2'], **args)
            if data["SMILES2"] is not None:
                if data["SMILES2"] !=  ".":
                    mod_compound.update(smiles=data["SMILES2"])
            if data["SMILES2"] is None:
                mod_compound.structure = None
            
        except Exception as e:
            raise e
            # if exception is of type value error, return the error message
            if type(e) == ValueError:
                return None, None, None, None, str(e)
            # if exception is of type other, return the error message
            else:
                return None, None, None, None, "Error loading compounds"
        
        if main_compound.structure is None:
            return None, None, None, None, "Error loading SMILES1"

        siteLocator = ModiFinder(main_compound, mod_compound, helpers=helper_compounds, **args)
        

        if mod_compound.structure is not None:
            if not (mod_compound.structure.HasSubstructMatch(main_compound.structure) or main_compound.structure.HasSubstructMatch(mod_compound.structure)):
                return None, None, None, None, "None of the structures are substructures of the other"
            if mod_compound.structure.HasSubstructMatch(main_compound.structure) and main_compound.structure.HasSubstructMatch(mod_compound.structure):
                return None, None, None, None, "Structures are the same"
        
        siriusText = "SIRIUS data was not available"
        # else:
        #     print("SIRIUS data was not available", data['USI1'])
        # if siteLocator.main_compound.Precursor_MZ > siteLocator.modified_compound.Precursor_MZ:
        #     return None, "Molecule precursor mass is higher than modified precursor mass", siriusText
        # else:
        args = copy.deepcopy(data)
        # remove SMILES and USI from args
        args.pop('SMILES1', None)
        args.pop('SMILES2', None)
        args.pop('USI1', None)
        args.pop('USI2', None)
        
        main_compound_peaks = [(main_compound.spectrum.mz[i], main_compound.spectrum.intensity[i]) for i in range(len(main_compound.spectrum.mz))]
        mod_compound_peaks = [(mod_compound.spectrum.mz[i], mod_compound.spectrum.intensity[i]) for i in range(len(mod_compound.spectrum.mz))]
        matched_peaks = siteLocator.get_edge_detail(main_compound.id, mod_compound.id)
        if matched_peaks is None:
            matched_peaks = []
        else:
            matched_peaks = matched_peaks.get_matches_pairs()
        peaksObj = {
            "main_compound_peaks": main_compound_peaks,
            "mod_compound_peaks": mod_compound_peaks,
            "matched_peaks": matched_peaks,
            "args": args,
            "main_precursor_mz": main_compound.spectrum.precursor_mz,
            "mod_precursor_mz": mod_compound.spectrum.precursor_mz,
        }

        fragmentsObj = {
            "frags_map": main_compound.spectrum.peak_fragments_map,
            "structure": main_compound.structure,
            "peaks": main_compound_peaks,
            "Precursor_MZ": main_compound.spectrum.precursor_mz,
        }

        return base64.b64encode(pickle.dumps(siteLocator)).decode(),  siriusText, base64.b64encode(pickle.dumps(peaksObj)).decode(), base64.b64encode(pickle.dumps(fragmentsObj)).decode(), None