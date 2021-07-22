"""
Description: Utilities for the 14N/15N calculator
"""

# Import packages
import math

import pandas as pd
import pyteomics
from pyopenms.pyopenms_5 import MzXMLFile, MSExperiment
from pyteomics import mass

modifications: dict = {}


def create_modification_list(mod_csv_filepath: str):
    """
    Add modifications to the dictionary.
    :param mod_csv_filepath: The path to the Modification CSV.
    :return:
    """
    # Read modification CSV
    mod_csv: pd.DataFrame = pd.read_csv(filepath_or_buffer=mod_csv_filepath)
    for _, row in mod_csv.iterrows():
        mass.std_aa_comp[f"{row['Prefix']}{row['Residue']}"] = mass.std_aa_comp[row['Residue']] + \
                                                               mass.Composition(formula=row['Composition'])
        modifications[(row['Residue'], row['Mass'])] = row['Prefix']


def create_modified_sequence(sequence: str, seq_start: int, mods: list) -> str:
    """
    Create the modified sequence.
    :param sequence: The sequence.
    :param seq_start: The position of the first residue to calculate the positions of modified residues.
    :param mods: The modifications with residue, position and mass change.
    :return: The modified sequence.
    """
    # Convert sequence to residues
    residues = list(sequence)
    for mod in mods:
        if (mod['type'], mod['modified']) in modifications:
            mod_index = mod['at'] - seq_start
            residues[mod_index] = f"{modifications[(mod['type'], mod['modified'])]}{residues[mod_index]}"
        else:
            print(f"\tWARN: {mod['modified']}@{mod['type']} is not supported.")

    # Join the residue to sequence. Hack: Convert ambiguous I/L (Letter: J) to I
    mod_sequence = "".join(residues).replace('J', 'I')
    return mod_sequence


def calculate_mass(seq_charge: tuple):
    """
    Calculate the 14N and 15N mass for a given sequence and the m/z for a given charge. :param seq_charge: Tuple
    containing the sequence and the charge. Also the constant minimum m/z-value to filter of too small peptides .
    :return: The tuple containing the 14N and 15N mass and m/z (for the given charge).
    """
    sequence: str = seq_charge[0]
    charge: int = seq_charge[1]
    minimum_mz: float = seq_charge[2]
    # Get the composition
    try:
        compositions = list(pyteomics.mass.isotopologues(sequence, elements_with_isotopes='N'))
        n14_comp = compositions[0]
        n15_comp = compositions[-1]
        # Calculate the mass
        n14_mass = round(pyteomics.mass.calculate_mass(n14_comp, type='M'), 3)
        n15_mass = round(pyteomics.mass.calculate_mass(n15_comp, type='M'), 3)
        # Get the m/z
        n14_mz = round(pyteomics.mass.calculate_mass(n14_comp, charge=charge), 3)
        n15_mz = round(pyteomics.mass.calculate_mass(n15_comp, charge=charge), 3)
        if n14_mz < minimum_mz or n15_mz < minimum_mz:
            return 0, 0, 0, 0
        else:
            return n14_mass, n14_mz, n15_mass, n15_mz
    except pyteomics.auxiliary.structures.PyteomicsError:
        print(f"\tCould not parse {seq_charge}")
        return 0, 0, 0, 0


def calculate_intensities(peptide_df: pd.DataFrame, mzxml_filepath: str, tolerance: float = 0.01) -> pd.DataFrame:
    """
    Calculate intensities from a mzXML file.
    :param peptide_df: The dataframe containing information about the peptide.
    :param mzxml_filepath: The path to the mzXML.
    :param tolerance: The tolerance in Da. Default is 0.01
    :return: The input dataframe with the intensity information.
    """
    # Load MS file
    exp = MSExperiment()
    MzXMLFile().load(mzxml_filepath, exp)
    # Remove all MSn spectra
    spec = []
    for s in exp.getSpectra():
        if s.getMSLevel() == 1:
            s.sortByPosition()
            spec.append(s)
    exp.setSpectra(spec)

    # Create lists to hold additional data
    scan_numbers: list = []
    n14_ints: list = []
    n15_ints: list = []
    n14_mzs: list = []
    n15_mzs: list = []
    ratios: list = []

    # Loop over the rows
    for _, row in peptide_df.iterrows():
        try:
            # Get spectrum from with the given RT
            spectrum = [s for s in exp if math.floor(row['RT']) < s.getRT() < math.ceil(row['RT'])][0]

            # Get the the closest peak closest to 14N and 15N m/z, respectively, within a given tolerance.
            n14_peak = spectrum[spectrum.findNearest(row['14N mz (Thr)'], tolerance / row['Charge'])]
            n15_peak = spectrum[spectrum.findNearest(row['15N mz (Thr)'], tolerance / row['Charge'])]
        except IndexError:
            print(f"\tWARN: Could not add {row['Sequence']}. No scans with N14 peak = {row['14N mz (Thr)']} m/z and "
                  f"and N15 peak = {row['14N mz (Thr)']} m/z at RT = {row['RT']} s.")
            scan_numbers.append(None)
            n14_mzs.append(None)
            n15_mzs.append(None)
            n14_ints.append(None)
            n15_ints.append(None)
            ratios.append(None)
            continue

        n14_int = n14_peak.getIntensity()
        n15_int = n15_peak.getIntensity()
        if n14_int == 0:
            ratio = None
        else:
            ratio = round(n15_int / n14_int, 3)

        # Set the data
        scan_numbers.append(spectrum.getNativeID().split('=')[1])
        n14_mzs.append(round(n14_peak.getMZ(), 3))
        n15_mzs.append(round(n15_peak.getMZ(), 3))
        n14_ints.append(round(n14_int, 3))
        n15_ints.append(round(n15_int, 3))
        ratios.append(ratio)

    # Add data
    peptide_df['Scan number'] = scan_numbers
    peptide_df['14N mz (Exp)'] = n14_mzs
    peptide_df['15N mz (Exp)'] = n15_mzs
    peptide_df['14N Int'] = n14_ints
    peptide_df['15N Int'] = n15_ints
    peptide_df['Ratio'] = ratios
    # Remove rows with no values
    peptide_df = peptide_df[peptide_df['Ratio'].notna()]

    return peptide_df
