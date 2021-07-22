"""
Description: Find hits in both 14N and 15N
"""

import pandas as pd
import math
import re
from tqdm import tqdm
from pyopenms.pyopenms_5 import MSExperiment, MzXMLFile


def read_peptide_lists(peptide_list_file: str, n14_tab_name: str, n15_tab_name: str) -> tuple:
    """
    Read the peptide list
    :param peptide_list_file: The filepath to the peptide list.
    :param n14_tab_name: The name of the 14N tab.
    :param n15_tab_name: The name of the 15N tab.
    :return: A tuple containing the 14N and 15N data.
    """
    used_columns = ['z', 'MH+ exp', 'MH+ theo', 'delta', 'from', 'to', 'seq', 'm/z', 'modifs']
    column_rename_dict: dict = {'z': 'Charge', 'MH+ exp': 'Mass (Exp)', 'MH+ theo': 'Mass (Thr)', 'delta': 'Delta',
                                'from': 'Start', 'to': 'End', 'seq': 'Sequence', 'modifs': 'Modifications'}
    new_column_order = ['Sequence', 'Charge', 'm/z', 'Mass (Exp)', 'Mass (Thr)', 'Start', 'End', 'Modifications']
    # Get, rename and order data
    n14_data: pd.DataFrame = pd.read_excel(peptide_list_file, sheet_name=n14_tab_name, usecols=used_columns,
                                           header=0)
    n14_data = n14_data.rename(columns=column_rename_dict)
    n14_data = n14_data[new_column_order]
    n15_data: pd.DataFrame = pd.read_excel(peptide_list_file, sheet_name=n15_tab_name, usecols=used_columns,
                                           header=0)
    n15_data = n15_data.rename(columns=column_rename_dict)
    n15_data = n15_data[new_column_order]
    return n14_data, n15_data


def find_hits(n14: pd.DataFrame, n15: pd.DataFrame, output_filepath: str):
    """
    Find the hits.
    :param n14: The 14N peptide list
    :param n15: The 15N peptide list
    :param output_filepath: The output filepath.
    :return: The hit list as a pandas dataframe
    """
    hit_list: pd.DataFrame = pd.DataFrame(columns=['Sequence', 'Modifications', 'Start', 'End', 'Charge',
                                                   '14N m/z', '15N m/z',
                                                   '14N Mass (Exp)', '15N Mass (Exp)'])
    # Loop each row in the peptide list
    for _, n14_row in tqdm(n14.iterrows(), total=n14.shape[0]):
        for _, n15_row in n15.iterrows():
            # Check if match
            if n14_row['Sequence'] == n15_row['Sequence'] and n14_row['Modifications'] == n15_row['Modifications'] and \
                    n14_row['Charge'] == n15_row['Charge']:
                # Add match to hit list
                hit_list = hit_list.append(pd.DataFrame([[n14_row['Sequence'], n14_row['Modifications'],
                                                          n14_row['Start'], n14_row['End'], n14_row['Charge'],
                                                          n14_row['m/z'], n15_row['m/z'], n14_row['Mass (Exp)'],
                                                          n15_row['Mass (Exp)']]],
                                                        columns=['Sequence', 'Modifications', 'Start', 'End',
                                                                 'Charge', '14N m/z', '15N m/z', '14N Mass (Exp)',
                                                                 '15N Mass (Exp)']),
                                           ignore_index=True)
    # Remove duplicates, sort sequences in alphabetical order and reset the index
    hit_list = hit_list.drop_duplicates(subset=['Sequence', 'Modifications'])
    hit_list = hit_list.sort_values(by=['Start'])
    hit_list = hit_list.reset_index(drop=True)
    hit_list.to_excel(output_filepath, sheet_name="List")

    return hit_list


if __name__ == '__main__':
    peptide_file = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\PeptideList.xlsx"
    hit_output_file = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\Hits\hits_37.xlsx"
    n14_tab = "Mix_37_14N"
    n15_tab = "Mix_37_15N"
    n14_list, n15_list = read_peptide_lists(peptide_list_file=peptide_file, n14_tab_name=n14_tab, n15_tab_name=n15_tab)
    hits: pd.DataFrame = find_hits(n14_list, n15_list, output_filepath=hit_output_file)
