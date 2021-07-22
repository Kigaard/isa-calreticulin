"""
Description: Find hits for between 14N and 15N
"""

# Import packages
import pandas as pd
from tqdm import tqdm


def read_peptide_lists(peptide_list_file: str, n14_tab_name: str, n15_tab_name: str) -> tuple:
    """
    Read the peptide list
    :param peptide_list_file: The filepath to the peptide list.
    :param n14_tab_name: The name of the 14N tab.
    :param n15_tab_name: The name of the 15N tab.
    :return: A tuple containing the 14N and 15N data.
    """
    # Get, rename and order data
    n14_data: pd.DataFrame = pd.read_excel(peptide_list_file, sheet_name=n14_tab_name, header=0)
    n14_data = n14_data[n14_data['V'] == 'Y']
    n15_data: pd.DataFrame = pd.read_excel(peptide_list_file, sheet_name=n15_tab_name, header=0)
    n15_data = n15_data[n15_data['V'] == 'Y']
    return n14_data, n15_data


def find_hits(n14: pd.DataFrame, n15: pd.DataFrame, output_filepath: str):
    """
        Find the hits.
        :param n14: The 14N peptide list
        :param n15: The 15N peptide list
        :param output_filepath: The output filepath.
        :return: The hit list as a pandas dataframe
        """
    hit_list: pd.DataFrame = pd.DataFrame(columns=n14.columns)
# Loop each row in the peptide list
    for _, n14_row in tqdm(n14.iterrows(), total=n14.shape[0]):
        for _, n15_row in n15.iterrows():
            # Check if match
            if n14_row['seq'] == n15_row['seq'] and n14_row['modifs'] == n15_row['modifs'] and \
                    n14_row['z'] == n15_row['z'] and n14_row['z'] == n15_row['z']:
                # Add match to hit list
                hit_list = hit_list.append(pd.DataFrame([n14_row], columns=n14.columns), ignore_index=True)

    hit_list = hit_list.drop_duplicates(subset=['seq', 'modifs'])
    hit_list = hit_list.sort_values(by=['from'])
    hit_list = hit_list.reset_index(drop=True)
    hit_list.to_excel(output_filepath, sheet_name="List")
    return hit_list


if __name__ == '__main__':
    peptide_file = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\PeptideList_new.xlsx"
    hit_output_file = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\Hits\hits_37.xlsx"
    n14_tab = "Mix_37_14N"
    n15_tab = "Mix_37_15N"
    n14_list, n15_list = read_peptide_lists(peptide_list_file=peptide_file, n14_tab_name=n14_tab, n15_tab_name=n15_tab)
    hits: pd.DataFrame = find_hits(n14_list, n15_list, output_filepath=hit_output_file)
