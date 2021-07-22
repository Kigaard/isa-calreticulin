"""
Description: Calculate the ratio between 14N and 15N peaks
"""
# Import packages
import time

from N145CalculatorUtilities import create_modified_sequence, create_modification_list, calculate_intensities, \
    calculate_mass
import pandas as pd
import pyteomics
from pyteomics import tandem


def read_tandem_result_file(tandem_result_filepath: str) -> pd.DataFrame:
    """
    Read the X!Tandem result file.
    :param tandem_result_filepath: The X!Tandem result file file path.
    :return: The peptide dataframe.
    """
    # Create dataframe
    peptide_df: pd.DataFrame = pd.DataFrame(columns=['Sequence', 'ModSequence', 'Charge', 'Start', 'End', 'RT',
                                                     '14N mass', '14N mz (Thr)', '15N mass', '15N mz (Thr)',
                                                     'Modifications'])
    # Read tandem result file
    tandem_result = pyteomics.tandem.filter(tandem_result_filepath, fdr=0.05)

    for result in tandem_result:
        for protein in result['protein']:
            # Handle peptide modification
            mod_sequence = protein['peptide']['seq']
            modification_str = None
            if 'aa' in protein['peptide']:
                mod_sequence = create_modified_sequence(protein['peptide']['seq'], protein['peptide']['start'],
                                                        protein['peptide']['aa'])

                modification_str = ", ".join([f"{mod['type']}{mod['at']}@{mod['modified']}" for mod
                                              in protein['peptide']['aa']])

            # Add row with sequence to the dataframe
            peptide_df = peptide_df.append(pd.DataFrame([[protein['peptide']['seq'], mod_sequence, result['z'],
                                                          protein['peptide']['start'], protein['peptide']['end'],
                                                          result['rt'], modification_str]],
                                                        columns=['Sequence', 'ModSequence', 'Charge', 'Start', 'End',
                                                                 'RT', 'Modifications']),
                                           ignore_index=True)
    peptide_df.to_excel('test.xlsx')

    return peptide_df


def setup_df(df: pd.DataFrame, min_mz: float) -> pd.DataFrame:
    # Clean the data
    df = df.drop_duplicates(subset=['ModSequence'], keep='first')
    df = df.reset_index(drop=True)
    df.to_excel('filter_test.xlsx')
    # Add region information
    df['StartEnd'] = df.apply(lambda x: (x['Start'], x['End']), axis=1)
    df['Region'] = df['StartEnd'].map(get_region)
    # Calculate the 14N and 15N masses along with m/z
    df['SeqCharge'] = df.apply(lambda x: (x['ModSequence'], x['Charge'], min_mz), axis=1)
    df['14N mass'], df['14N mz (Thr)'], df['15N mass'], df['15N mz (Thr)'] = \
        zip(*df['SeqCharge'].map(calculate_mass))
    df = df[df['15N mz (Thr)'] != 0]

    # Remove unnecessary columns
    del df['SeqCharge']
    del df['ModSequence']
    del df['StartEnd']

    return df


def get_region(start_end: tuple) -> str:
    start_region = get_region_from_position(start_end[0])
    end_region = get_region_from_position(start_end[1])
    if start_region == end_region:
        return start_region
    else:
        return f"{start_region};{end_region}"


def get_region_from_position(pos: int) -> str:
    if pos < 17:
        # No peptide should be in this region
        return "Signal"
    elif pos < 204:
        return "Core"
    elif pos < 305:
        return "P"
    elif pos < 336:
        return "Core"
    else:
        return "C"


# Entry point
if __name__ == '__main__':
    tolerance = 0.01
    minimum_mz = 500

    tandem_file = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\XTandem\EXP3_01258_VM_mix_37.xml"
    mzxml_file = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\mzXML\EXP3_01258_VM.mzXML"
    excel_result_file = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\ratios.xlsx"

    start_time = time.time()

    # Create modification list
    print("Create modification list")
    create_modification_list("Modifications.csv")

    # Read the X!Tandem result file
    print("Read X!Tandem file")
    peptides_raw: pd.DataFrame = read_tandem_result_file(tandem_file)

    # Setup dataframe
    print("Setup dataframe")
    peptides: pd.DataFrame = setup_df(peptides_raw, minimum_mz)

    # Calculate the intensities
    print("Calculate intensities")
    peptides = calculate_intensities(peptides, mzxml_file, tolerance=tolerance)

    # Reorder columns
    print("Reorder columns before creating Excel-file")
    peptides = peptides[['Sequence', 'Charge', 'Start', 'End', 'Modifications', 'RT', 'Scan number', '14N mass',
                         '14N mz (Thr)', '14N mz (Exp)', '14N Int', '15N mass', '15N mz (Thr)', '15N mz (Exp)',
                         '15N Int', 'Ratio', 'Region']]
    peptides = peptides.reset_index(drop=True)

    # Write the data to the Excel-file
    print("Write result to Excel-file")
    peptides.to_excel(excel_result_file)

    run_time = time.time() - start_time
    print("Done in {} seconds ({} minutes).".format(round(run_time, 2), round(run_time/60, 2)))

