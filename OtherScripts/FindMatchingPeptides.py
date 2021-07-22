"""
Description: Find peptide which matches across the different conditions
"""
# Import packages
import pandas as pd
from tqdm import tqdm


def find_matching_hits(hits37: pd.DataFrame, hits42: pd.DataFrame, hits42zn: pd.DataFrame, output: str) -> pd.DataFrame:
    hit_list: pd.DataFrame = pd.DataFrame(columns=['Sequence', 'Modifications', 'Start', 'End', 'Charge',
                                                   '14N m/z', '15N m/z',
                                                   '14N Mass (Exp)', '15N Mass (Exp)'])
    for _, hit37 in tqdm(hits37.iterrows(), total=hits37.shape[0]):
        for _, hit42 in hits42.iterrows():
            if hit37['Sequence'] == hit42['Sequence'] and \
                    hit37['Modifications'] == hit42['Modifications'] and \
                    hit37['Start'] == hit42['Start']:
                for _, hit42zn in hits42zn.iterrows():
                    if hit37['Sequence'] == hit42['Sequence'] == hit42zn['Sequence'] and \
                            hit37['Modifications'] == hit42['Modifications'] == hit42zn['Modifications'] and \
                            hit37['Start'] == hit42['Start'] == hit42zn['Start']:
                        hit_list = hit_list.append(pd.DataFrame([[hit37['Sequence'], hit37['Modifications'],
                                                                  hit37['Start'], hit37['End'], hit37['Charge'],
                                                                  hit37['14N m/z'], hit37['15N m/z'],
                                                                  hit37['14N Mass (Exp)'], hit37['15N Mass (Exp)']]],
                                                                columns=['Sequence', 'Modifications', 'Start', 'End',
                                                                         'Charge', '14N m/z', '15N m/z',
                                                                         '14N Mass (Exp)',
                                                                         '15N Mass (Exp)']),
                                                   ignore_index=True)
                hit_list = hit_list.drop_duplicates(subset=['Sequence', 'Modifications'])
                hit_list = hit_list.sort_values(by=['Start'])
                hit_list = hit_list.reset_index(drop=True)
                hit_list.to_excel(output, sheet_name="List")


if __name__ == '__main__':
    hits37_df: pd.DataFrame = pd.read_excel(io=r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\Hits\hits_37.xlsx",
                                            sheet_name="List")
    hits42_df: pd.DataFrame = pd.read_excel(io=r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\Hits\hits_42.xlsx",
                                            sheet_name="List")
    hits42_zn_df: pd.DataFrame = pd.read_excel(io=r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\Hits\hits_42_Zn"
                                                  r".xlsx",
                                               sheet_name="List")
    output_file: str = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\Hits\MatchingHits.xlsx"
    matching_hits: pd.DataFrame = find_matching_hits(hits37=hits37_df, hits42=hits42_df, hits42zn=hits42_zn_df,
                                                     output=output_file)
