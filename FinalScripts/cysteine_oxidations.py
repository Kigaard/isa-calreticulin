"""
Description: Find the cysteine oxidation and plot them in a bar chart.
"""

# Import packages
import os.path

import pandas as pd

from quantiative_plot_utilities import create_plots_from_peptide_lists

"""
The modification dictionary and the cysteine positions, which are specific for CRT.
"""
MODIFICATIONS: dict = {
    -33.988: "Dehydro",
    15.995: "Oxidation",
    31.990: "Sulfinic",
    63.962: "SulfDiOx",
    47.967: "SulfOx",
    47.985: "Sulfonic",
    79.957: "SSulfonic",
    91.957: "SO3"
}
POSITIONS = [105, 137, 163]


def _combine_and_clean_modifications(raw_df: pd.DataFrame) -> pd.DataFrame:
    """
    Combine 'SulfOx' and 'Sulfinic acid' as they have the same mass and perform some cleaning up after the
    combination.

    :param raw_df: The  raw data frame.
    :return: The combined and cleaned data frame
    """
    raw_df['SulfinicSulfOx'] = raw_df['Sulfinic'] + raw_df['SulfOx']
    del raw_df['Sulfinic']
    del raw_df['SulfOx']
    return raw_df[['Dehydro', 'Oxidation', 'Sulfonic', 'SulfDiOx', 'SulfinicSulfOx', 'SSulfonic', 'SO3']]


if __name__ == '__main__':
    BASE_FILE_PATH = r""
    conditions_n14 = [['Tryp_rCrt14_Cys', 'A) Trypsin rCrt14 (rCRT)'], ['Tryp_pCrt14_Cys', 'B) Trypsin pCrt14 (pCRT)'],
                      ['Mix_37', 'C) 37 °C (pCRT)'], ['Mix_42', 'D) 42 °C (pCRT)'],
                      ['Mix_42_Zn', 'E) 42 °C + Zn (pCRT)']]
    conditions_n15 = [['Tryp_rCrt14_Cys_N15', 'A) Trypsin rCrt14 (rCRT)'],
                      ['Tryp_pCrt14_Cys_N15', 'B) Trypsin pCrt14 (rCRT)'],
                      ['Mix_37_N15', 'C) 37 °C (rCRT)'], ['Mix_42_N15', 'D) 42 °C (rCRT)'],
                      ['Mix_42_Zn_N15', 'E) 42 °C + Zn (rCRT)']]

    labels: list = ["Dehydroalanine (-34)", "Oxidation (16)", "Sulfinic acid (32)",
                    "Persulfinic acid (64)", "Sulfonic acid / SulfOx (48)", "Persulfonic acid (80)", "SO3 (92)"]
    
    # Create the plot for 14N data
    file_data_n14: list = []
    for condition in conditions_n14:
        file_data_n14.append((os.path.join(BASE_FILE_PATH, condition[0]), condition[1], None))
    create_plots_from_peptide_lists(peptide_lists=file_data_n14, modifications=MODIFICATIONS,
                                    modification_position=POSITIONS,
                                    combine_function=_combine_and_clean_modifications, labels=None)
    
    # Create the plot for 15N data
    file_data_n15: list = []
    for condition in conditions_n15:
        file_data_n15.append((os.path.join(BASE_FILE_PATH, condition[0]), condition[1], None))
    create_plots_from_peptide_lists(peptide_lists=file_data_n15, modifications=MODIFICATIONS,
                                    modification_position=POSITIONS,
                                    combine_function=_combine_and_clean_modifications, labels=None)
