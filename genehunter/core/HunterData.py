from collections import OrderedDict

import os
import sys
import time
import shutil
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import h5py as h5
import pandas as pd
import numpy as np

from genehunter.core.mtcorr import get_bh_thres


class GwasData(object):
    def __init__(self):
        self._data = None

    @property
    def data(self):
        return self._data

    # @data.setter
    # def data(self, df):
    #     self._data = df

    def read_hdf5(self, filepath, pval_threshold=1.0, mac_threshold=0, fdr_alpha=0.05):
        file_df = None
        sys.stdout.write("Reading file: {} ... ".format(os.path.basename(filepath)))
        sys.stdout.flush()
        score_threshold = -np.log10(pval_threshold)
        with h5.File(filepath, "r") as h5file:
            start = time.time()
            root = h5file["pvalues"]

            all_pvalues = []
            all_selected = []
            groupnames = root.keys()
            for gname in groupnames:
                grp = root[gname]

                raw_pvalues = np.power(10.0, -grp["scores"].value)
                all_pvalues.extend(raw_pvalues)
                # if self._raw_pvalues is not None:
                #     self._raw_pvalues.append(raw_pvalues)
                # else:
                #     self._raw_pvalues = raw_pvalues

                mth = grp["macs"][...] >= mac_threshold
                pth = grp["scores"][...] >= score_threshold
                combinedth = mth & pth
                all_selected.extend(combinedth)

                nrow = combinedth.sum()
                if nrow > 0:
                    # tmpdf = pd.DataFrame({
                    #     "origin": np.repeat(os.path.basename(filepath), nrow),
                    #     "chromosomes": np.repeat(int(gname.strip("chr")), nrow),
                    #     "positions": grp["positions"][combinedth],
                    #     "pvalues": raw_pvalues[combinedth],
                    #     "macs": grp["macs"][combinedth]
                    # })
                    tmpdf = pd.DataFrame(columns=["Original_file", "Chromosome", "SNP_pos", "GWAS_pvalue", "MAC"])
                                         #dtype={'origin': np.str, "chromosomes": np.int32, "positions": np.int64, "pvalues": np.float64, "macs": np.int32})
                    tmpdf["Original_file"] = np.repeat(os.path.basename(filepath), nrow).astype(np.str)
                    # tmpdf["chromosomes"] = np.repeat(int(gname.strip("chr")), nrow).astype(np.int32)
                    tmpdf["Chromosome"] = np.repeat(gname, nrow)
                    tmpdf["SNP_pos"] = grp["positions"][combinedth].astype(np.int64)
                    # pvals = raw_pvalues[combinedth].tolist()
                    tmpdf["GWAS_pvalue"] = raw_pvalues[combinedth]
                    tmpdf["MAC"] = grp["macs"][combinedth].astype(np.int32)

                    if file_df is not None:
                        file_df = pd.concat([file_df, tmpdf], ignore_index=True, axis=0)
                    else:
                        file_df = tmpdf
        sys.stdout.write("[ {:f}s ]\n".format(time.time() - start))

        if file_df.shape[0] > 0:
            sys.stdout.write("{:d} positions passed.\n".format(file_df.shape[0]))
            sys.stdout.write("calculating thresholds ... ")
            sys.stdout.flush()
            start = time.time()
            bonferroni_th = fdr_alpha / len(all_pvalues)
            bh_rejected, bh_corrected = fdrcorrection0(all_pvalues, fdr_alpha, method="indep")
            bh_th = get_bh_thres(all_pvalues, fdr_alpha)
            file_df["Bonferroni_{:.3f}_threshold".format(fdr_alpha)] = np.repeat(bonferroni_th, file_df.shape[0])
            file_df["BH_{:.3f}_threshold".format(fdr_alpha)] = np.repeat(bh_th["thes_pval"], file_df.shape[0])
            file_df["BH_FDR_{:.3f}_adjusted".format(fdr_alpha)] = np.array(bh_corrected)[all_selected]
            file_df["BH_FDR_{:.3f}_rejected".format(fdr_alpha)] = np.array(bh_rejected)[all_selected]
            sys.stdout.write("[ {:f}s ]\n\n".format(time.time() - start))
            if self._data is not None:
                pd.concat([self._data, file_df], ignore_index=True, axis=0)
            else:
                self._data = file_df
        else:
            sys.stdout.write("no positions passed.")


    @staticmethod
    def remap_hdf5_results(h5filepath, mapfilepath):
        chrom_mapping = dict()
        with open(mapfilepath, "r") as mapfile:
            for line in mapfile:
                cols = line.split(',')
                try:
                    chrom_mapping[int(cols[0])] = cols[1].strip()
                except ValueError:
                    continue

        with h5.File(h5filepath, "a") as h5file:
            pval_grp = h5file["pvalues"]

            backup_filepath = h5filepath + ".bak"
            if not os.path.exists(backup_filepath):
                shutil.copy2(h5filepath, backup_filepath)
            else:
                sys.stdout.write("Backup copy for file: {} exists. Skipping this file.\n".format(h5filepath))
                return

            current_keys = pval_grp.keys()
            if len(set(chrom_mapping.values()) & set(current_keys)) == len(current_keys):
                sys.stdout.write("Chromosomes seem to be already remapped in file: {}. Skipping this file.\n".format(h5filepath))
                return

            sys.stdout.write("Remapping chromosomes in file: {}.\n".format(h5filepath))
            for key in pval_grp.keys():
                intkey = int(key.strip('chr'))
                real_name = chrom_mapping[intkey]
                pval_grp.move(key, real_name)


class InputData(object):
    def __init__(self, fdr=0.05):
        self.fdr = fdr
        self.df = pd.DataFrame(columns=["Original_file", "Chromosome", "SNP_pos", "GWAS_p-value", "uDist", "dDist"])
                                               #FDR_{}_rejected".format(fdr), "FDR_{}_adjusted_p-value".format(fdr),
                                               #"Bonferroni_{}_threshold".format(fdr), "BH_{}_threshold".format(fdr),
                                               #"BHY_{}_threshold".format(fdr), "Gene_start", "Gene_end",
                                               #"Gene_orientation", "Relative_Distance", "SNP_relative_position",
                                               #"target_AGI", "target_element_type", "target_sequence_type",
                                               #"target_annotation", "target_attributes"])
        self.cols = self.df.columns

    def add_dataset(self, data):
        self.df = pd.concat([self.df, data], axis=0, ignore_index=True)
        self.df = self.df.reindex_axis(self.cols, axis=1)
        # self.df.append(data, ignore_index=True)


class OutputData(object):
    def __init__(self, fdr=0.05):
        self.fdr = fdr
        self.col_types = OrderedDict([
            ("Original_file", "str"),
            ("Chromosome", "str"),
            ("SNP_pos", "int"),
            ("GWAS_p-value", "float"),
            ("FDR_{}_rejected".format(fdr), "float"),
            ("FDR_{}_adjusted_p-value".format(fdr), "float"),
            ("Bonferroni_{}_threshold".format(fdr), "float"),
            ("BH_{}_threshold".format(fdr), "float"),
            ("BHY_{}_threshold".format(fdr), "float"),
            ("Gene_start", "int"),
            ("Gene_end", "int"),
            ("Gene_orientation", 'S1'),
            ("Relative_Distance", "int"),
            ("SNP_relative_position", "str"),
            ("target_AGI", "str"),
            ("target_element_type", "str"),
            ("target_sequence_type", "str"),
            ("target_annotation", "str"),
            ("target_attributes", "str")
            ])
        self.df = pd.DataFrame({k: pd.Series(dtype=v) for k, v in self.col_types.items()})

    def add_dataset(self, dataseries):
        self.df = self.df.append(dataseries, ignore_index=True)

    def write_csv(self, filepath):
        self.df = self.df.reindex_axis(self.col_types.keys(), axis=1)
        self.df.to_csv(filepath, index=False, na_rep="NA")


if __name__ == '__main__':
    gwd = GwasData()
    gwd.read_hdf5("/net/gmi.oeaw.ac.at/busch/lab/Marco/GWAS/Medicago/20170606_Mt_StantonGeddes2013_gwas-results/20170606_Mt_StantonGeddes2013_height3.hdf5", pval_threshold=1.0e-6, mac_threshold=10)
    gwd.read_hdf5("/net/gmi.oeaw.ac.at/busch/lab/Marco/GWAS/Medicago/20170606_Mt_StantonGeddes2013_gwas-results/20170606_Mt_StantonGeddes2013_totNod.hdf5", pval_threshold=1.0e-6, mac_threshold=10)