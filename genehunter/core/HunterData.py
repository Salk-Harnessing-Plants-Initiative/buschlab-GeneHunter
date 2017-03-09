from collections import OrderedDict

import pandas as pd
import numpy as np

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
