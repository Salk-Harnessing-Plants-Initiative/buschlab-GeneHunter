from collections import OrderedDict

import os
import sys
import time
import h5py as h5
import pandas as pd
import numpy as np


class GwasData(object):
    def __init__(self):
        self._data = None

    @property
    def data(self):
        return self._data

    # @data.setter
    # def data(self, df):
    #     self._data = df

    def read_hdf5(self, filepath, pval_threshold=1.0, mac_threshold=0):
        sys.stdout.write("Reading file: {} ... ".format(os.path.basename(filepath)))
        sys.stdout.flush()
        score_threshold = -np.log10(pval_threshold)
        with h5.File(filepath, "r") as h5file:
            start = time.time()
            root = h5file["pvalues"]

            all_pvalues = []
            groupnames = root.keys()
            for gname in groupnames:
                grp = root[gname]

                raw_pvalues = np.power(10.0, -grp["scores"].value)
                all_pvalues.append(raw_pvalues)
                # if self._raw_pvalues is not None:
                #     self._raw_pvalues.append(raw_pvalues)
                # else:
                #     self._raw_pvalues = raw_pvalues

                mth = grp["macs"][...] >= mac_threshold
                pth = grp["scores"][...] >= score_threshold
                combinedth = mth & pth

                nrow = combinedth.sum()
                if nrow > 0:
                    # tmpdf = pd.DataFrame({
                    #     "origin": np.repeat(os.path.basename(filepath), nrow),
                    #     "chromosomes": np.repeat(int(gname.strip("chr")), nrow),
                    #     "positions": grp["positions"][combinedth],
                    #     "pvalues": raw_pvalues[combinedth],
                    #     "macs": grp["macs"][combinedth]
                    # })
                    tmpdf = pd.DataFrame(columns=["origin", "chromosomes", "positions", "pvalues", "macs"])
                                         #dtype={'origin': np.str, "chromosomes": np.int32, "positions": np.int64, "pvalues": np.float64, "macs": np.int32})
                    tmpdf["origin"] = np.repeat(os.path.basename(filepath), nrow).astype(np.str)
                    tmpdf["chromosomes"] = np.repeat(int(gname.strip("chr")), nrow).astype(np.int32)
                    tmpdf["positions"] = grp["positions"][combinedth].astype(np.int64)
                    # pvals = raw_pvalues[combinedth].tolist()
                    tmpdf["pvalues"] = raw_pvalues[combinedth]
                    tmpdf["macs"] = grp["macs"][combinedth].astype(np.int32)

                    if self._data is not None:
                        self._data = pd.concat([self._data, tmpdf], ignore_index=True, axis=0)
                    else:
                        self._data = tmpdf
        sys.stdout.write("[ {:f}s ]\n".format(time.time() - start))
        sys.stdout.write("calculating thresholds ... ")
        sys.stdout.flush()


        print(self._data)



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