import sys
import traceback
import requests
import json
from typing import List, Tuple

from modifinder import Compound

adduct_mapping = {'M+H': '[M+H]+',
'[M+H]': '[M+H]+',
'[M+H]+': '[M+H]+',
'[M+H]1+': '[M+H]+',
'M+H]': '[M+H]+',
'M+Na': '[M+Na]+',
'[M+Na]': '[M+Na]+',
'[M+Na]+': '[M+Na]+',
'[M+Na]1+': '[M+Na]+',
'2M+Na': '[2M+Na]+',
'M2+Na': '[2M+Na]+',
'[2M+Na]+': '[2M+Na]+',
'[2M+Na]': '[2M+Na]+',
'M+K': '[M+K]+',
'[M+K]': '[M+K]+',
'[M+K]+': '[M+K]+',
'[M+K]1+': '[M+K]+',
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
'[M-H]1-': '[M-H]-',
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
'[M+NH4]1+': '[M+NH3+H]+',
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
'[M+Cl]1-': '[M+Cl]-',
'M+K-2H': '[M-2H+K]-',
'[M-2H+K]': '[M-2H+K]-',
'M-2H]': '[M-2H]2-',
'M-2H': '[M-2H]2-',
'M-2H-': '[M-2H]2-',
'M+Na-2H': '[M-2H+Na]-',
'[M-2H+Na]': '[M-2H+Na]-',
'M+Br': '[M+Br]-',
'[M+Br]1-': '[M+Br]-',
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

def filter_peaks_by_ratio_to_base_peak(peaks:List[Tuple[float,float]], ratio_to_base_peak:float = 0.01):
        """Remove peaks with intensity lower than a given ratio to the base peak.
        
        Parameters
        ----------
        peaks : List[Tuple[float, float]]
            List of (mz, intensity) tuples representing the spectrum peaks.
        ratio_to_base_peak : float (0, 1), default is 0.01
            The ratio to the base peak.
        """

        base_peak_intensity = max([intensity for (mz, intensity) in peaks])
        new_mz = []
        new_intensity = []
        for index, intensity in enumerate([peak[1] for peak in peaks]):
            if intensity >= float(ratio_to_base_peak) * base_peak_intensity:
                new_mz.append(peaks[index][0])
                new_intensity.append(intensity)
        
        return list(zip(new_mz, new_intensity))
    


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
            data = _get_from_metabolomics_resolver(identifier)
            data['id'] = identifier
            data = _harmonize_spectrum_keys(data)

            # Sort peaks if needed
            if 'peaks' in data and isinstance(data['peaks'], list) and len(data['peaks']) > 0:
                data['peaks'] = sorted(data['peaks'], key=lambda x: x[0])
            return data

    link = "https://external.gnps2.org/gnpsspectrum?SpectrumID={}".format(identifier)
    try:
        res = requests.get(link)
        parsed = res.json()
    except Exception:
        data = _get_from_metabolomics_resolver(original_identifier)
        data['usi'] = original_identifier
        data['id'] = identifier
        data = _harmonize_spectrum_keys(data)
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

    data = _harmonize_spectrum_keys(data)
    data['id'] = identifier

    # Ensure peaks are sorted
    if 'peaks' in data and isinstance(data['peaks'], list) and len(data['peaks']) > 0:
        data['peaks'] = sorted(data['peaks'], key=lambda x: x[0])
    return data

def load_helpers(
        data: List[str],
        ratio_to_base_peak: float = None,
        ) -> List[Compound]:
    """ Load helpers from a list of identifiers, failing gracefully if the smile string is invalid.
    """
    loaded_helpers = []
    failed_helpers = []
    for h in data:
        try:
            lh = get_data(h)

            if ratio_to_base_peak:
                lh['peaks'] = filter_peaks_by_ratio_to_base_peak(lh['peaks'], ratio_to_base_peak=ratio_to_base_peak)

            ch = Compound(
                spectrum=lh['peaks'],
                precursor_mz=lh['precursor_mz'],
                precursor_charge=lh['precursor_charge'],
                adduct=lh.get('adduct', None),
                smiles=lh.get('smiles', None)
            )

            loaded_helpers.append(ch)
        except Exception as e:
            # Print the traceback
            print(f"Error loading helper compound {h}: {str(e)}", flush=True)
            traceback.print_exc(file=sys.stderr)
            failed_helpers.append(h)
            raise e
    
    print(f"Loaded {len(loaded_helpers)} helper compounds successfully. Failed to load {len(failed_helpers)} helper compounds: {failed_helpers}", flush=True)
    return loaded_helpers

def _harmonize_spectrum_keys(data):
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

def _get_from_metabolomics_resolver(identifier: str) -> dict:
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

    data = _harmonize_spectrum_keys(data)
    return data