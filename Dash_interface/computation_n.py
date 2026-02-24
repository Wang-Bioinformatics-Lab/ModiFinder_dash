import traceback

from dash import Dash, html, dcc, Input, Output, State, dash_table
import base64
import pickle
import json
import copy
from urllib.parse import quote
from typing import List
import sys

import requests
from modifinder import ModiFinder, Compound
from rdkit import Chem
from datetime import datetime

adduct_mapping = {'M+H': '[M+H]+',
'[M+H]': '[M+H]+',
'[M+H]+': '[M+H]+',
'M+H]': '[M+H]+',
'M+Na': '[M+Na]+',
'[M+Na]': '[M+Na]+',
'[M+Na]+': '[M+Na]+',
'2M+Na': '[2M+Na]+',
'M2+Na': '[2M+Na]+',
'[2M+Na]+': '[2M+Na]+',
'[2M+Na]': '[2M+Na]+',
'M+K': '[M+K]+',
'[M+K]': '[M+K]+',
'[M+K]+': '[M+K]+',
'[2M+K]+': '[2M+K]+',
'2M+K': '[2M+K]+',
'[2M+K]': '[2M+K]+',
'M+H-H20': '[M-H2O+H]+',
'M+H-H2O': '[M-H2O+H]+',
'[M-H2O+H]+': '[M-H2O+H]+',
'M-H20+H': '[M-H2O+H]+',
'[M+H-H2O]+': '[M-H2O+H]+',
'M-H2O+H': '[M-H2O+H]+',
'M+H-2H2O': '[M-2H2O+H]+',
'M-2H2O+H': '[M-2H2O+H]+',
'[M-2H2O+H]+': '[M-2H2O+H]+',
'M-2(H2O)+H': '[M-2H2O+H]+',
'2M+Na-2H': '[2M-2H+Na]-',
'2M-2H+Na': '[2M-2H+Na]-',
'M-H': '[M-H]-',
'[M-H]': '[M-H]-',
'[M-H]-': '[M-H]-',
'M-H-': '[M-H]-',
'M-H1': '[M-H]-',
'3M+Na': '[3M+Na]+',
'[3M+Na]+': '[3M+Na]+',
'[M]+': '[M]+',
'M+': '[M]+',
'M-e': '[M]+',
'M2+H': '[2M+H]+',
'2M+H': '[2M+H]+',
'[2M+H]+': '[2M+H]+',
'[2M+H]': '[2M+H]+',
'[M+2H]': '[M+2H]2+',
'[M+2H]2+': '[M+2H]2+',
'M+2H]': '[M+2H]2+',
'M+2H+2': '[M+2H]2+',
'M+2H': '[M+2H]2+',
'M+acetate': '[M+CH3COOH-H]-',
'M+CH3COOH-H': '[M+CH3COOH-H]-',
'M+CH3COO': '[M+CH3COOH-H]-',
'M+ACN+H': '[M+CH3CN+H]+',
'[M+ACN+H]+': '[M+CH3CN+H]+',
'[M+H+CH3CN]': '[M+CH3CN+H]+',
'M+2Na': '[M+2Na]2+',
'M+2Na]': '[M+2Na]2+',
'M+HCOO': '[M+HCOOH-H]-',
'[M-H+HCOOH]': '[M+HCOOH-H]-',
'M+FA-H': '[M+HCOOH-H]-',
'M+formate': '[M+HCOOH-H]-',
'[M+H+HCOOH]': '[M+HCOOH-H]-',
'2M+FA-H': '[2M+HCOOH-H]-',
'[2M-H+HCOOH]': '[2M+HCOOH-H]-',
'M+NH4': '[M+NH3+H]+',
'[M+NH4]+': '[M+NH3+H]+',
'[M+NH4]': '[M+NH3+H]+',
'2M+Hac-H': '[2M+CH3COOH-H]-',
'2M-H': '[2M-H]-',
'[2M-H]': '[2M-H]-',
'2M+NH4': '[2M+NH3+H]+',
'[2M+NH4]+': '[2M+NH3+H]+',
'[2M+NH4]': '[2M+NH3+H]+',
'[2M+Ca]2+': '[2M+Ca]2+',
'[M+Ca]2+': '[M+Ca]2+',
'[3M+Ca]2+': '[3M+Ca]2+',
'[2M+Ca-H]+': '[2M-H+Ca]+',
'[2M-H2O+H]+': '[2M-H2O+H]+',
'[4M+Ca]2+': '[4M+Ca]2+',
'[3M+NH4]+': '[3M+NH3+H]+',
'3M+NH4': '[3M+NH3+H]+',
'[2M-2H2O+H]+': '[2M-2H2O+H]+',
'[M+ACN+NH4]+': '[M+CH3CN+NH3+H]+',
'[5M+Ca]2+': '[5M+Ca]2+',
'[3M+K]+': '[3M+K]+',
'[3M+Ca-H]+': '[3M-H+Ca]2+',
'[M-H+2Na]+': '[M-H+2Na]+',
'M-H+2Na': '[M-H+2Na]+',
'[M-3H2O+H]+': '[M-3H2O+H]+',
'M-3H2O+H': '[M-3H2O+H]+',
'[M-3H2O+2H]2+': '[M-3H2O+2H]2+',
'[M-2H2O+2H]2+': '[M-2H2O+2H]2+',
'[M-4H2O+H]+': '[M-4H2O+H]+',
'[M-5H2O+H]+': '[M-5H2O+H]+',
'[M+Ca-H]+': '[M+Ca-H]+',
'[2M-H+2Na]+': '[2M-H+2Na]+',
'[2M-3H2O+H]+': '[2M-3H2O+H]+',
'[M+H+Na]2+': '[M+Na+H]2+',
'[M-2H2O+NH4]+': '[M-2H2O+NH3+H]+',
'[2M-2H+Na]': '[2M-2H+Na]-',
'[M-H+CH3OH]': '[M+CH3OH-H]-',
'M+MeOH-H': '[M+CH3OH-H]-',
'M-H2O-H': '[M-H2O-H]-',
'[M-H-H2O]': '[M-H2O-H]-',
'M+Cl-': '[M+Cl]-',
'M+Cl': '[M+Cl]-',
'[M+Cl]': '[M+Cl]-',
'M+K-2H': '[M-2H+K]-',
'[M-2H+K]': '[M-2H+K]-',
'M-2H]': '[M-2H]2-',
'M-2H': '[M-2H]2-',
'M-2H-': '[M-2H]2-',
'M+Na-2H': '[M-2H+Na]-',
'[M-2H+Na]': '[M-2H+Na]-',
'M+Br': '[M+Br]-',
'3M-H': '[3M-H]-',
'[3M-H]': '[3M-H]-',
'[M+H+CH3OH]': '[M+CH3OH+H]+',
'M+CH3OH+H': '[M+CH3OH+H]+',
'[2M+H+CH3CN]': '[2M+CH3CN+H]+',
'M-CO2-H': '[M-CO2-H]-',
'[2M-2H+K]': '[2M-2H+K]-',
'2M+K-2H': '[2M-2H+K]-',
'[M+Na+CH3CN]': '[M+CH3CN+Na]+',
'M-H2+H': '[M-H2+H]-',
'M-H+Cl]': '[M-H+Cl]2-',
'M-H+Cl': '[M-H+Cl]2-',
'3M+H': '[3M+H]+',
'[3M+H]': '[3M+H]+',
'M+H-NH3': '[M-NH3+H]+',
'M-NH3+H': '[M-NH3+H]+',
'M-H+C2H2O': '[M+C2H2O-H]-',
'M+H-C2H2O': '[M+C2H2O-H]-',
'M-H+CH2O2': '[M+CH2O2-H]-',
'M+CH2O2-H': '[M+CH2O2-H]-',
'M+TFA-H': '[M+C2HF3O2-H]-',
'M-C2HF3O2-H': '[M+C2HF3O2-H]-',
'[M]1+': '[M]1+'}


gnps_keys_mapping = {
    ## precursor
    "precursor_mz": "precursor_mz",
    ## charge
    "precursor_charge": "precursor_charge",
    "charge": "precursor_charge", 
    ## smiles
    "smiles": "smiles",
    "smile": "smiles",
    ## adduct
    "adduct": "adduct",
    ## peaks
    "peaks": "peaks",
    ## instrument
    "instrument": "instrument",
    ## name
    "name": "name",
    "compound_name": "name",
    ## spectrum_id
    "spectrum_id": "spectrum_id",
    "spectrumid": "spectrum_id",
    ## exact mass
    "exact_mass": "exact_mass",
    "exactmass": "exact_mass",
    ## mz
    "fragment_mz": "mz",
    "mz": "mz",
    "mzs": "mz",
    ## intensity
    "fragment_intensities": "intensity",
    "intensities": "intensity",
}

def filter_peaks_by_ratio_to_base_peak(spectrum, ratio_to_base_peak:float = 0.01):
        """Remove peaks with intensity lower than a given ratio to the base peak.
        
        Parameters
        ----------
        ratio_to_base_peak : float (0, 1), default is 0.01
            The ratio to the base peak.
        change_spectrum : bool, default is True
            If True, the peaks with intensity lower than the given ratio will be removed in place.
            If False, a new Spectrum object with the peaks removed will be returned.
        """
        
        base_peak = max(spectrum.intensity)
        new_mz = []
        new_intensity = []
        for index, intensity in enumerate(spectrum.intensity):
            if intensity >= float(ratio_to_base_peak) * base_peak:
                new_mz.append(spectrum.mz[index])
                new_intensity.append(intensity)
        
        spectrum.mz = new_mz
        spectrum.intensity = new_intensity

        return spectrum

def remove_larger_than_precursor_peaks(spectrum):
        """
        Remove peaks that are larger than the precursor m/z value.
        """
        
        new_mz = []
        new_intensity = []
        for mz, intensity in zip(spectrum.mz, spectrum.intensity):
            if mz < spectrum.precursor_mz * 0.99:
                new_mz.append(mz)
                new_intensity.append(intensity)
        
        spectrum.mz = new_mz
        spectrum.intensity = new_intensity

        return spectrum


    
def harmonize_spectrum_keys(data):
    """
    Parse the data to a universal format.

    This function takes a dictionary of data and converts it into a universal format.
    It processes specific keys like "peaks_json" and "Charge" differently, and attempts
    to convert other values to floats. If the conversion to float is successful and the
    key is "Charge", it further converts the value to an integer.

    Args:
        :data (dict): The input data dictionary to be parsed.

    Returns:
        :dict: A dictionary with keys converted to a universal format and values processed
              accordingly.
    """
    def _convert_to_universal_key(key: str) -> str:
        """
        Convert different types of keys to universal keys.
        This function standardizes various key names to a universal format. 

        Args:
            :key (str): The key to be converted.
        
        Returns:
            :str: The converted key.
        """
        key = key.lower()
        key = key.replace(" ", "_")
        return gnps_keys_mapping.get(key, key)

    res = {}
    for key, value in data.items():
        converted_key = _convert_to_universal_key(key)
        if key == "peaks_json":
            res['peaks'] = json.loads(value)
        elif converted_key == "adduct":
            res[converted_key] = adduct_mapping.get(value, value)
        else:
            try:
                if converted_key in ["precursor_charge", "precursor_mz", "ms_level", "scan", "exact_mass"]:
                    value = float(value)
                if converted_key in ["precursor_charge", "charge", "ms_level"]:
                        value = int(value)
            except Exception:
                raise ValueError(f"Could not convert {key} to number")
            res[converted_key] = value
    return res

def get_from_metabolomics_resolver(identifier: str) -> dict:
    """
    Get partial data (ms2 data) from USI
    param identifier: str - USI
    return: dict - dictionary of data with keys: precursor_mz, precursor_charge, mz: list, intensity: list
    """
    url = 'https://metabolomics-usi.gnps2.org/json/' + "?usi1=" + identifier
    try:
        r = requests.get(url)
        data = json.loads(r.text)
    except:
        raise Exception("Error in retrieving data from GNPS for identifier: {}, link: {}".format(identifier, url))

    data = harmonize_spectrum_keys(data)
    return data

def get_data(identifier: str) -> dict:
    """
    Get data from GNPS, either from USI or Accession. if the identifier points to a known item in gnps,
      it will return the full data, otherwise it will return partial data (ms2 data)
    param identifier: str - USI or Accession
    return: dict - dictionary of data
    """

    data = dict()
    data['usi'] = None

    if "mzspec" in identifier:                              # It's a USI
        data['usi'] = identifier

        if "accession" in identifier:                       #       It's a library spectrum
            original_identifier = str(identifier)
            identifier = identifier.split(":")[-1]
        else:                                               #       It's a USI that isn't a library spectrum
            data = get_from_metabolomics_resolver(identifier)
            data['id'] = identifier
            data = harmonize_spectrum_keys(data)

            # Sort peaks if needed
            if 'peaks' in data and isinstance(data['peaks'], list) and len(data['peaks']) > 0:
                data['peaks'] = sorted(data['peaks'], key=lambda x: x[0])

            return data

    link = "https://external.gnps2.org/gnpsspectrum?SpectrumID={}".format(identifier)
    try:
        res = requests.get(link)
        parsed = res.json()
    except Exception:
        data = get_from_metabolomics_resolver(original_identifier)
        data['usi'] = original_identifier
        data['id'] = identifier
        data = harmonize_spectrum_keys(data)
        return data

    try:
        data.update(parsed['annotations'][0])
    except KeyError:
        pass
    try:
        data.update(parsed['spectruminfo'])
    except KeyError:
        pass
    try:
        data['comments'] = parsed['comments']
    except KeyError:
        pass

    data = harmonize_spectrum_keys(data)
    data['id'] = identifier

    # Ensure peaks are sorted
    if 'peaks' in data and isinstance(data['peaks'], list) and len(data['peaks']) > 0:
        data['peaks'] = sorted(data['peaks'], key=lambda x: x[0])

    return data

def load_helpers(
        data: List[str],
        ratio_to_base_peak: float = None,
        remove_large_peaks: bool = True
        ) -> List[Compound]:
    """ Load helpers from a list of identifiers, failing gracefully if the smile string is invalid.
    """
    loaded_helpers = []
    failed_helpers = []
    for h in data:
        try:
            lh = get_data(h)
            ch = Compound(
                spectrum=lh['peaks'],
                precursor_mz=lh['precursor_mz'],
                precursor_charge=lh['precursor_charge'],
                adduct=lh.get('adduct', None),
                smiles=lh.get('smiles', None)
            )
            if ratio_to_base_peak:
                ch.spectrum = filter_peaks_by_ratio_to_base_peak(ch.spectrum, ratio_to_base_peak=ratio_to_base_peak)
            if remove_large_peaks:
                ch.spectrum = remove_larger_than_precursor_peaks(ch.spectrum)
            loaded_helpers.append(ch)
        except Exception as e:
            # Print the traceback
            print(f"Error loading helper compound {h}: {str(e)}", flush=True)
            traceback.print_exc(file=sys.stderr)
            failed_helpers.append(h)
            raise e
    
    print(f"Loaded {len(loaded_helpers)} helper compounds successfully. Failed to load {len(failed_helpers)} helper compounds: {failed_helpers}", flush=True)
    return loaded_helpers

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
        usi1 = args.pop('USI1', None)
        usi2 = args.pop('USI2', None)
        
        spectrum1 = get_data(usi1)
        spectrum2 = get_data(usi2)
        if spectrum1['adduct'] is None:
            # Replace with adduct from data
            spectrum1['adduct'] = data.get('adduct', None)
        if spectrum2['adduct'] is None:
            # Replace with adduct from data
            spectrum2['adduct'] = data.get('adduct', None)

        # TODO: What to do if adduct differs at this point?

        # TODO: Filter adducts in Helpers?

        # Options propagated out of ModiFinder 
        ratio_to_base_peak = args.pop('filter_peaks_variable', None)
        remove_large_peaks = True
        
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
            remove_large_peaks=remove_large_peaks
        )

        if data["SMILES1"] == "" or data["SMILES1"] is None:
            data["SMILES1"] = None
        
        if data["SMILES2"] == "" or data["SMILES2"] is None:
            data["SMILES2"] = None

        try:
            if data['adduct']:
                args['adduct'] = data['adduct']
            main_compound = Compound(
                spectrum=spectrum1['peaks'],
                precursor_mz=spectrum1['precursor_mz'],
                precursor_charge=spectrum1['precursor_charge'],
                adduct=spectrum1['adduct'],
                smiles=data["SMILES1"]
            )
            mod_compound = Compound(
                spectrum=spectrum2['peaks'],
                precursor_mz=spectrum2['precursor_mz'],
                precursor_charge=spectrum2['precursor_charge'],
                adduct=spectrum2['adduct'],
                smiles=data["SMILES2"] if data["SMILES2"] is not None and data["SMILES2"] != "" else None
            )
            
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
        
        # Perform actions for  ratio_to_base_peak filter
        if ratio_to_base_peak:
            ratio_to_base_peak = float(ratio_to_base_peak)
            main_compound.spectrum = filter_peaks_by_ratio_to_base_peak(main_compound.spectrum, ratio_to_base_peak)
            mod_compound.spectrum = filter_peaks_by_ratio_to_base_peak(mod_compound.spectrum, ratio_to_base_peak)

        # Perform actions for remove_large_peaks filter
        if remove_large_peaks:
            main_compound.spectrum = remove_larger_than_precursor_peaks(main_compound.spectrum)
            mod_compound.spectrum = remove_larger_than_precursor_peaks(mod_compound.spectrum)

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