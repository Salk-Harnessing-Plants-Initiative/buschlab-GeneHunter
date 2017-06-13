from collections import OrderedDict

import h5py as h5
import pandas as pd
import numpy as np


class GwasData(object):
    def __init__(self):
        self.df = None;

    def read_hdf5(self, filepath, pval_threshold=1.0, mac_threshold=0):
        self.df = None

        with h5.File(filepath, "r") as h5file:
            root = h5file["pvalues"]

            groupnames = root.keys()
            for gname in groupnames:
                tmpdf = pd.DataFrame()
                grp = root[gname]
                mth = grp["macs"][...] >= mac_threshold
                pth = grp["scores"][...] >= -np.log10(pval_threshold)

                combinedth = mth & pth

                tmpdf["positions"] = grp["positions"][combinedth]
                tmpdf["pvalues"] = grp["scores"][combinedth]
                tmpdf["macs"] = grp["macs"][combinedth]
                tmpdf.insert(0, "chromosomes", np.repeat(int(gname.strip("chr")), tmpdf.shape[0]))

                if self.df is not None:
                    pd.concat([self.df, tmpdf], ignore_index=True)
                else:
                    self.df = tmpdf
        pass




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
    gwd.read_hdf5("/net/gmi.oeaw.ac.at/busch/lab/Marco/GWAS/Medicago/20170606_Mt_StantonGeddes2013_gwas-results/20170606_Mt_StantonGeddes2013_height3.hdf5", pval_threshold=1.0e-5, mac_threshold=0)