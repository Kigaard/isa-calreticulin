"""
Description: Find the methionine and proline and plot them in a bar chart.
"""

# Import packages
import os.path

from quantiative_plot_utilities import create_plots_from_peptide_lists

"""
The modification dictionary and the cysteine positions, which are specific for CRT.
"""
MODIFICATIONS: dict = {
    15.995: "Oxidation"
}
POSITIONS = [19, 83, 114, 122, 131, 134, 139, 178, 204, 205, 211, 216, 228, 233, 239, 243, 245, 250, 257, 263, 264, 269,
             277, 283, 292, 297, 301, 303, 357, 410]

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
                                    combine_function=None, labels=None)
    
    # Create the plot for 15N data
    file_data_n15: list = []
    for condition in conditions_n15:
        file_data_n15.append((os.path.join(BASE_FILE_PATH, condition[0]), condition[1], None))

    create_plots_from_peptide_lists(peptide_lists=file_data_n15, modifications=MODIFICATIONS,
                                    modification_position=POSITIONS,
                                    combine_function=None, labels=None)
