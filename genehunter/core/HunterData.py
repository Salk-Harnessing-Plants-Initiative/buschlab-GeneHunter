import pandas as pd

class InputData(object):
    def __init__(self, fdr):
        self.fdr = fdr
        self.inputdata = pd.DataFrame(columns=["Original_file", "Chromosome", "SNP_pos", "GWAS_p-value", "uDist", "dDist"])
                                               #FDR_{}_rejected".format(fdr), "FDR_{}_adjusted_p-value".format(fdr),
                                               #"Bonferroni_{}_threshold".format(fdr), "BH_{}_threshold".format(fdr),
                                               #"BHY_{}_threshold".format(fdr), "Gene_start", "Gene_end",
                                               #"Gene_orientation", "Relative_Distance", "SNP_relative_position",
                                               #"target_AGI", "target_element_type", "target_sequence_type",
                                               #"target_annotation", "target_attributes"])

    def add_dataset(self, dataseries):
        self.inputdata.append(dataseries, ignore_index=True)


class OutputData(object):
    def __init__(self, fdr):
        self.inputdata = pd.DataFrame(columns=["Original_file", "Chromosome", "SNP_pos", "GWAS_p-value",
                                               "FDR_{}_rejected".format(fdr), "FDR_{}_adjusted_p-value".format(fdr),
                                               "Bonferroni_{}_threshold".format(fdr), "BH_{}_threshold".format(fdr),
                                               "BHY_{}_threshold".format(fdr), "Gene_start", "Gene_end",
                                               "Gene_orientation", "Relative_Distance", "SNP_relative_position",
                                               "target_AGI", "target_element_type", "target_sequence_type",
                                               "target_annotation", "target_attributes"])
