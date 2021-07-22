"""
Description: Calculate hit intensities from 14N/15N hits
"""

# Import packages
import pandas as pd
from tqdm import tqdm
from pyopenms.pyopenms_5 import MSExperiment, MzXMLFile, MSSpectrum


def calculate_intensity(experiment: MSExperiment, n14_mz: float, n15_mz,
                        sequence: str, charge: int, modification: str, mod_seq: str) -> pd.DataFrame:
    hits: pd.DataFrame = pd.DataFrame(columns=['Sequence', 'Modifications', 'Charge', 'RT', 'Scan number',
                                               '14N m/z (Exp)', '14N m/z (Thr)', '14N Intensity',
                                               '14N m/z (Exp)', '15N m/z (Thr)', '15N Intensity',
                                               'Ratio', 'ModSeq'])

    tolerance = 0.01
    for scan in experiment:
        # Get peaks
        n14_peak_index = scan.findNearest(n14_mz, tolerance / charge)
        n15_peak_index = scan.findNearest(n15_mz, tolerance / charge)
        if n14_peak_index != -1 and n15_peak_index != -1:
            n14_peak = scan[n14_peak_index]
            n15_peak = scan[n15_peak_index]

            n14_int = n14_peak.getIntensity()
            n15_int = n15_peak.getIntensity()

            if n14_int == 0 and n15_int == 0:
                continue

            if n14_int == 0:
                ratio = 0
            else:
                ratio = round(n15_int / n14_int, 3)

            hits = hits.append(
                pd.DataFrame([[sequence, modification, charge, scan.getRT(), scan.getNativeID().split('=')[1],
                               n14_mz, n14_peak.getMZ(), n14_int, n15_mz, n15_peak.getMZ(), n15_int, ratio, mod_seq]],
                             columns=hits.columns), ignore_index=True)

    return hits


def read_ms_data(mzxml_filepath: str):
    # Read mzXML file
    print("Read mzXML file...")
    exp: MSExperiment = MSExperiment()
    MzXMLFile().load(mzxml_filepath, exp)
    spec = []
    for s in exp.getSpectra():
        if s.getMSLevel() == 1:
            spec.append(s)
    exp.setSpectra(spec)
    print("Done reading")
    return exp


def calculate_intensities(hit_list_filepath: str, mzxml_filepath: str,
                          intensity_hit_list_filepath: str) -> pd.DataFrame:
    hit_list: pd.DataFrame = pd.read_excel(hit_list_filepath, header=0)
    ms_experiment = read_ms_data(mzxml_filepath=mzxml_filepath)

    intensity_hitlist: pd.DataFrame = pd.DataFrame(columns=['Sequence', 'Modifications', 'Charge', 'RT', 'Scan number',
                                                            '14N m/z (Exp)', '14N m/z (Thr)', '14N Intensity',
                                                            '14N m/z (Exp)', '15N m/z (Thr)', '15N Intensity',
                                                            'Ratio', 'ModSeq'])
    hit_list['ModSeq'] = hit_list.apply(lambda x: create_mod_sequence_string(x['Sequence'], x['Modifications'],
                                                                             x['Start']), axis=1)

    for _, hit in tqdm(hit_list.iterrows(), total=hit_list.shape[0]):
        hit_scans: pd.DataFrame = calculate_intensity(experiment=ms_experiment, n14_mz=float(hit['14N m/z']),
                                                      n15_mz=float(hit['15N m/z']), sequence=hit['Sequence'],
                                                      charge=int(hit['Charge']), modification=hit['Modifications'],
                                                      mod_seq=hit['ModSeq'])
        intensity_hitlist = intensity_hitlist.append(hit_scans, ignore_index=True)

    intensity_hitlist = intensity_hitlist.set_index(keys=['ModSeq'], append=False)
    intensity_hitlist.to_excel(intensity_hit_list_filepath)
    return intensity_hitlist


def create_mod_sequence_string(sequence: str, modifications: str, peptide_start: int) -> str:
    if modifications == "-":
        return sequence

    mod_sequence = ""
    for mod in modifications.split(" "):
        resid = int(mod.split('@')[0])
        mass = float(mod.split('@')[1])
        relative_resid = resid - peptide_start
        mod_sequence = sequence[:relative_resid+1] + f"({mass})" + sequence[relative_resid+1:]

    return mod_sequence


if __name__ == '__main__':
    hit_list_path = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\Hits\hits_trypsin_rCRT14.xlsx"
    mzxml_path = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\mzXML\EXP3_01353_VM_tryp_mix_rCrt14.mzXML"
    intensity_hit_list_path = r"C:\Users\Mads\Desktop\ISA_Spring_2021\Autodigestion\HitsIntensity" \
                              r"\Hit_intensity_Trypsin_rCrt14.xlsx"

    hit_df = calculate_intensities(hit_list_filepath=hit_list_path, mzxml_filepath=mzxml_path,
                                   intensity_hit_list_filepath=intensity_hit_list_path)
