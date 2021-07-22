"""
Description: Utility functions for creating the quantitative plots.
"""

import math
from typing import List, Tuple, Dict, Callable, Union

import pandas as pd
import matplotlib.pyplot as plt


def _find_modifications(hits_df: pd.DataFrame, positions: List[int], modification_dict: Dict[float, str]) \
        -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Find the modifications in the given peptide DataFrame.

    :param hits_df: The hits found hits.
    :param positions: The list of positions available for modification.
    :param modification_dict: The dictionary with the modification mass and the name.
    :return: The tuple containing DataFrames with position, modifications for the absolute and percentage values,
        and the total peptide count, respectively.
    """
    # Get the modifications  and filter them
    modifications: list = _parse_and_filter_modifications(modifications=hits_df['modifs'], positions=positions,
                                                          modification_dict=modification_dict)
    # Count the number of modifications at each position.
    modification_count: dict = {}
    for pos in positions:
        modification_count[str(pos)]: dict = {}
        for modification_type in modification_dict.values():
            modification_count[str(pos)][modification_type]: int = 0

    for pos, mod in modifications:
        modification_count[pos][mod] += 1

    # Calculate the modification percentage
    modification_percentage: dict = {}

    peptide_count: dict = _count_peptides(hits_df=hits_df, positions=positions)

    for pos in modification_count:
        modification_percentage[pos]: dict = {}
        for mod in modification_count[pos]:
            if peptide_count[int(pos)] != 0:
                modification_percentage[pos][mod] = (modification_count[pos][mod] / peptide_count[int(pos)]) * 100
            else:
                modification_percentage[pos][mod] = 0

    # Convert the dictionaries to data frames
    percentage_df: pd.DataFrame = pd.DataFrame.from_dict(modification_percentage, orient='index')
    count_df: pd.DataFrame = pd.DataFrame.from_dict(modification_count, orient='index')
    peptide_count_df: pd.DataFrame = pd.DataFrame.from_dict(data=peptide_count, orient='index')
    return count_df, percentage_df, peptide_count_df


def _parse_and_filter_modifications(modifications: pd.Series, positions: List[int],
                                    modification_dict: Dict[float, str]) -> List[tuple]:
    """
    Parse and filter the modifications.

    :param modifications: The Series containing the modifications found in the hits.
    :param positions: The list of positions.
    :return: The list containing tuples with the modification position and the modification type.
    """
    mods_raw = [mod for mod in list(modifications) if mod != '-']
    mod_temp: list = []
    # Get all the modifications even if the multiple exist for one peptide
    for mod_raw in mods_raw:
        if ' ' in mod_raw:
            mods = str.split(' ')
            for mod in mods:
                mod_temp.append(mod)
        else:
            mod_temp.append(mod_raw)

    # Split the modifications into residue nad mass
    modifications_list = [str.split(mod, sep='@') for mod in mod_temp]
    # Get the modification position and type
    return [(mod[0], modification_dict[float(mod[1])]) for mod in modifications_list if int(mod[0]) in positions]


def _count_peptides(hits_df: pd.DataFrame, positions: List[int]) -> dict:
    """
    Count the the number of peptide that contains each of the specified cysteines.

    :param hits_df: The data frame containing the hits.
    :param positions: The list of positions available for modification.
    :return: The number of peptide for each position.
    """
    counts: dict = {}

    for pos in positions:
        counts[pos] = 0

    for start_pos, end_pos in zip(hits_df['from'], hits_df['to']):
        for pos in positions:
            if start_pos <= int(pos) <= end_pos:
                counts[pos] += 1

    return counts


def _create_plot(mod_dict: Dict[str, Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]], labels: Union[List[str], None],
                 max_y: int):
    """
    Create plot.

    :param mod_dict: The dictionary with the condition name and the count and percentage, and total count DataFrames.
    :param labels: The list of labels. If None, no legend will be shown.
    :param max_y: The maximum y-value.
    """
    # Create the subplots
    fig, ax = plt.subplots(nrows=len(mod_dict), sharex='all', figsize=(20 * 1 / 2.54, 20 * 1 / 2.54))
    fig.text(0.03, 0.025, 'Position', ha='center', va='center')  # X-label
    fig.text(0.02, 0.5, 'Percentage (Count)', ha='center', va='center', rotation='vertical')  # Y-label
    # Loop over the conditions and create a plot for each condition
    for condition_index, condition in enumerate(mod_dict):
        # Get the DataFrames
        count_df: pd.DataFrame = mod_dict[condition][0]
        percentage_df: pd.DataFrame = mod_dict[condition][1]
        peptide_count_df: pd.DataFrame = mod_dict[condition][2]

        # Create the bar plot
        percentage_df.plot.bar(ax=ax[condition_index], legend=False)
        ax[condition_index].set_title(condition)
        ax[condition_index].set_ylim([0, max_y])
        # Add annotations for each patch
        for patch_idx, p in enumerate(ax[condition_index].patches):
            # Get dataframe index
            row_idx = patch_idx % len(count_df.index)
            column_idx = math.floor(patch_idx / len(count_df.index))
            # Get the values
            percentage = percentage_df.iloc[row_idx, column_idx]
            count = count_df.iloc[row_idx, column_idx]
            total_count = peptide_count_df.iloc[row_idx, 0]
            label_y_offset: int = 20 if patch_idx % 2 == 0 else 0
            label: str = f"{format(percentage, '.1f').rstrip('0').rstrip('.')}%\n({count}/{total_count})" \
                if percentage != 0 else ''
            ax[condition_index].annotate(label, (p.get_x(), p.get_y()), xytext=(0, 10 + label_y_offset),
                                         textcoords='offset points')

    # Possibly create legend
    if labels is not None:
        fig.legend(loc='lower center', ncol=len(labels), labels=labels)
    plt.tight_layout(pad=0.5)
    plt.show()


def create_plots_from_peptide_lists(peptide_lists: List[Tuple[str, str, str]], modifications: Dict[float, str],
                                    modification_position: List[int], combine_function: Union[Callable, None],
                                    labels: Union[List[str], None], max_y: int = 100):
    """
    Create plots from the a list of peptide lists

    :param peptide_lists: The list of tuples containing the name of the peptide list and the condition and the
        sheet name. If sheet name is None, 'Sheet1' is used.
    :param modifications: The dictionary with the modification mass and the modification name.
    :param modification_position: The list of position available for modification.
    :param combine_function: The function which can be used for combining columns etc.
    :param labels: The labels to be used in the plot.
    :param max_y: The maximum y-value shown in the plot. Default 100.
    """
    modification_files: Dict[str, Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {}
    for peptide_list in peptide_lists:
        hits: pd.DataFrame = pd.read_excel(io=f"{peptide_list[0]}.xlsx",
                                           sheet_name=peptide_list[2] if peptide_list[2] is not None else 'Sheet1',
                                           usecols=["V", "modifs", "from", "to", "seq"])

        # Remove invalid peptides
        hits = hits[hits['V'] == "Y"]
        # Get the mods and the modifications percentages.
        mod_df_raw, percentage_df_raw, peptide_count_df = _find_modifications(hits_df=hits,
                                                                              positions=modification_position,
                                                                              modification_dict=modifications)

        # Combine if a combine function is given.
        if combine_function is not None:
            mod_df = combine_function(mod_df_raw)
            percentage_df = combine_function(percentage_df_raw)
        else:
            mod_df = mod_df_raw
            percentage_df = percentage_df_raw

        mod_df.index = [int(idx) for idx in mod_df.index]
        mod_df = mod_df.sort_index()
        percentage_df.index = [int(idx) for idx in percentage_df.index]
        percentage_df = percentage_df.sort_index()
        peptide_count_df.index = [int(idx) for idx in peptide_count_df.index]
        peptide_count_df = peptide_count_df.sort_index()

        modification_files[peptide_list[1]] = (mod_df, percentage_df, peptide_count_df)

    # Create the plots
    _create_plot(mod_dict=modification_files, labels=labels, max_y=max_y)
