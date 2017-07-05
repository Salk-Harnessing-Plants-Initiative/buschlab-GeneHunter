# -*- coding: utf-8 -*-

"""genehunter.genehunter: module doing the argument parsing and dispatching the commands."""

import sys
import os
import re
import glob
import argparse
import pandas as pd
import numpy as np
import statsmodels.sandbox.stats.multicomp as sm
import genehunter.core.mtcorr as mt
from genehunter.core.GeneAnnotationDbCreator import GeneAnnotationDbCreator
from genehunter.core.GeneAnnotationDbExtractor import GeneAnnotationDbExtractor
from genehunter.core.HunterData import InputData, GwasData

__version__ = "0.2.0"


def main():
    sys.stdout.write("gene annotation database - version %s\n" % __version__)
    GeneAnnotator()


class GeneAnnotator(object):
    def __init__(self):
        self.parse_arguments()
        self.inputdata = None


    @staticmethod
    def create_db(args):
        if os.path.exists(args.output) and not args.add:
            raise Exception("database file '%s' exists!\nChoose another output file. Terminating!\n" % args.output)
        if args.gff is not None and not os.path.exists(args.gff):
            raise Exception("cannot open GFF file '%s'! Terminating!\n" % args.gff)
        if args.desc is not None and not os.path.exists(args.desc):
            raise Exception("cannot open descriptions file '%s'! Terminating!\n" % args.desc)
        if args.aliases is not None and not os.path.exists(args.aliases):
            raise Exception("cannot open aliases file '%s'! Terminating!\n" % args.aliases)
        if args.sorf is not None and not os.path.exists(args.sorf):
            raise Exception("cannot open sORF file '%s'! Terminating!\n" % args.sorf)
        creator = GeneAnnotationDbCreator(args.output)
        creator.create_gff_entries(args.gff)
        # if args.desc is not None:
        #     creator.create_desc_entries(args.desc)
        # if args.aliases is not None:
        #     creator.create_name_entries(args.aliases)
        # if args.sorf is not None:
        #     creator.create_sorf_entries(args.sorf)
        return

    @staticmethod
    def print_stats(args):
        if not os.path.exists(args.db):
            raise Exception("cannot open database '%s'!" % args.db)

        extractor = GeneAnnotationDbExtractor(args.db)
        extractor.print_stats()

    def extract_loc(self, args):
        intervals = []
        self.inputdata = InputData()
        if not os.path.exists(args.db):
            raise Exception("cannot open database '%s'!" % args.db)

        if args.file is not None:
            in_df = pd.read_csv(args.file, sep=',', header=None)
            # with open(args.file, 'r') as ifile:
            #     for line in ifile:
            #         dseries = pd.Series()
            #         cols = re.split(',|\t', line)
            #         chrom = cols[0].strip()
            #         pos = int(cols[1].strip())
            #         dseries["Chromosome"] = chrom
            #         dseries["SNP_pos"] = pos
            #         if len(cols) == 2:
            #             dseries["uDist"] = args.loc1
            #             dseries["dDist"] = args.loc2
            #             # intervals.append((chrom, pos - args.loc1, pos + args.loc2))
            #         else:
            #             dseries["uDist"] = int(cols[2].strip())
            #             dseries["dDist"] = int(cols[3].strip())
            #             # intervals.append((chrom, pos - int(cols[2].strip()), pos + int(cols[3].strip())))
            self.inputdata.df["Chromosome"] = in_df.iloc[:, 0].astype(np.str)
            self.inputdata.df["SNP_pos"] = in_df.iloc[:, 1].astype(np.int32)
            if in_df.shape[1] > 2:
                self.inputdata.df["uDist"] = in_df.iloc[:, 2].astype(np.int32)
                self.inputdata.df["dDist"] = in_df.iloc[:, 3].astype(np.int32)
            else:
                self.inputdata.df["uDist"] = args.loc1
                self.inputdata.df["dDist"] = args.loc2
        elif args.c:
            dseries = pd.Series()
            dseries["Chromosome"] = args.chr
            dseries["SNP_pos"] = args.loc1
            dseries["uDist"] = args.loc2
            dseries["dDist"] = args.loc2
            self.inputdata.add_dataset(dseries)
            # intervals.append((args.chr, args.loc1 - args.loc2, args.loc1 + args.loc2))
        elif args.i:
            dseries = pd.Series()
            dseries["Chromosome"] = args.chr
            dseries["SNP_pos"] = (args.loc1 + args.loc2) / 2
            dseries["uDist"] = (args.loc2 - args.loc1) / 2
            dseries["dDist"] = (args.loc2 - args.loc1) / 2
            self.inputdata.add_dataset(dseries)
            # intervals.append((args.chr, args.loc1, args.loc2))

        extractor = GeneAnnotationDbExtractor(args.db)
        outdf = extractor.extract_by_loc(self.inputdata)
        outdf.write_csv(args.output)

        # for interval in intervals:
        #     extractor.extract_by_loc(interval[1], interval[2], interval[0])
        # extractor.write_results(args.output, args.depth)

        return

    @staticmethod
    def extract_agi(args):
        if args.agi is None and args.file is None:
            sys.stdout.write("no AGI(s) given. Terminating!\n")
            exit(1)

        agis = []
        if args.file is not None:
            with open(args.file, 'r') as ifile:
                for line in ifile:
                    agis.append(line.strip())
        else:
            agis.append(args.agi)

        extractor = GeneAnnotationDbExtractor(args.db)
        for agi in agis:
            extractor.extract_by_agi(agi)
        extractor.write_results(args.output, 1)
        return

    @staticmethod
    def extract_hunter(args):
        dbextract = GeneAnnotationDbExtractor(args.db)
        sys.stdout.write("gene hunter using database: {}".format(args.db))

        all_peaks_df = None
        for gwasfilename in glob.glob(os.path.join(args.dir, args.name)):
            assert gwasfilename.endswith(".hdf5") #TODO: create stable test
            gwd = GwasData()
            gwd.read_hdf5(gwasfilename, pval_threshold=args.pvalue_threshold, mac_threshold=args.minor_allele_count, fdr_alpha=args.fdr)
            # peaks_df = gwd.data

            # genes_df = none
            for ix, row in gwd.data.iterrows():
                gw_chr = str(row["Chromosome"]).strip()
                gw_pos = row["SNP_pos"]
                genes = dbextract.extract_loc_uddist(gw_chr, gw_pos, args.udistance, args.ddistance)
                sys.stdout.write("    peak: {}, pos {} -> {} genes in range\n".format(gw_chr, gw_pos, len(genes)))
                if len(genes) == 0:
                    ext_row = pd.Series(row)
                    ext_row["Gene_start"] = "NA"
                    ext_row["Gene_end"] = "NA"
                    ext_row["Gene_orientation"] = "NA"
                    ext_row["Relative_distance"] = "NA"
                    ext_row["SNP_relative_position"] = "NA"
                    ext_row["Target_AGI"] = "NA"
                    ext_row["Target_element_type"] = "NA"
                    ext_row["Target_sequence_type"] = "NA"
                    ext_row["Target_attributes"] = "NA"

                    if all_peaks_df is not None:
                        all_peaks_df = pd.concat([all_peaks_df, ext_row.to_frame().transpose()], axis=0, ignore_index=True)
                    else:
                        all_peaks_df = ext_row.to_frame().transpose()
                    continue

                for g in genes:
                    ext_row = pd.Series(row)
                    ext_row["Gene_start"] = g.start
                    ext_row["Gene_end"] = g.end
                    ext_row["Gene_orientation"] = g.strand
                    if g.strand == '+':
                        ext_row["Relative_distance"] = gw_pos - g.start
                    else:
                        ext_row["Relative_distance"] = g.start - gw_pos

                    if g.start <= gw_pos <= g.end:
                        ext_row["SNP_relative_position"] = "in gene"
                    elif gw_pos < g.start:
                        if g.strand == '+':
                            ext_row["SNP_relative_position"] = "upstream"
                        else:
                            ext_row["SNP_relative_position"] = "downstream"
                    else:
                        if g.strand == '+':
                            ext_row["SNP_relative_position"] = "downstream"
                        else:
                            ext_row["SNP_relative_position"] = "upstream"
                    ext_row["Target_AGI"] = g.id
                    ext_row["Target_element_type"] = g.feature
                    ext_row["Target_sequence_type"] = g.sequencetype
                    ext_row["Target_attributes"] = g.attribute

                    if all_peaks_df is not None:
                        all_peaks_df = pd.concat([all_peaks_df, ext_row.to_frame().transpose()], axis=0, ignore_index=True)
                    else:
                        all_peaks_df = ext_row.to_frame().transpose()

                    if args.depth >= 1:
                        for rna in g.rna:
                            ext_row = pd.Series(row)
                            ext_row["Gene_start"] = rna.start
                            ext_row["Gene_end"] = rna.end
                            ext_row["Gene_orientation"] = rna.strand
                            if rna.strand == '+':
                                ext_row["Relative_distance"] = gw_pos - rna.start
                            else:
                                ext_row["Relative_distance"] = rna.start - gw_pos

                            if rna.start <= gw_pos <= rna.end:
                                ext_row["SNP_relative_position"] = "in feature"
                            elif gw_pos < rna.start:
                                if rna.strand == '+':
                                    ext_row["SNP_relative_position"] = "upstream"
                                else:
                                    ext_row["SNP_relative_position"] = "downstream"
                            else:
                                if rna.strand == '+':
                                    ext_row["SNP_relative_position"] = "downstream"
                                else:
                                    ext_row["SNP_relative_position"] = "upstream"
                            ext_row["Target_AGI"] = rna.id
                            ext_row["Target_element_type"] = rna.feature
                            ext_row["Target_sequence_type"] = rna.sequencetype
                            ext_row["Target_attributes"] = rna.attribute

                            all_peaks_df = pd.concat([all_peaks_df, ext_row.to_frame().transpose()], axis=0, ignore_index=True)
            sys.stdout.write("\n")
        if args.output is not None:
            args.output = args.output.replace("_", "-")
            out_path = "{}_gene-hunter_u{:d}_d{:d}_pval{:.3e}_mac{:d}_fdr{:.3f}.txt".format(args.output, args.udistance, args.ddistance,
                                                                                            args.pvalue_threshold,
                                                                                            args.minor_allele_count,
                                                                                            args.fdr)
            out_path = os.path.join(args.dir, out_path)
            all_peaks_df.to_csv(out_path, sep='\t', header=True, index=False)
        else:
            all_peaks_df.to_string(sys.stdout, header=True, index=False)

    @staticmethod
    def extract_hunter_deprecated(args):
        dbextract = GeneAnnotationDbExtractor(args.db)
        dbid = dbextract.get_database_id()

        if dbid == 'Lj'.encode():
            sys.stdout.write("using the Lotus annotation database.\n")
        elif dbid == 'AT'.encode():
            sys.stdout.write("using the Arabidopsis Tair database.\n")
        else:
            sys.stdout.write("using unkown gene database (ID = {}). ".format(dbid))

        if args.output is None:
            ostream = sys.stdout
            isfile = False
        else:
            ostream = open(args.output, 'w')
            isfile = True

        ostream.write("Original_file\tChromosome\tSNP_pos\tGWAS_p-value\tFDR_{}_rejected\t".format(args.fdr))
        ostream.write("FDR_{}_adjusted_p-value\tBonferroni_{}_threshold\t".format(args.fdr, args.fdr))
        ostream.write("BH_{}_threshold\tBHY_{}_threshold\tGene_start\tGene_end\t".format(args.fdr, args.fdr, args.fdr))
        ostream.write("Gene_orientation\tRelative_Distance\tSNP_relative_position\ttarget_AGI\ttarget_element_type\t")
        ostream.write("target_sequence_type\ttarget_annotation\ttarget_attributes\n")
        # ostream.write("short_symbol\tlong_symbol\tshort_description\tlong_description\n")

        for gwasfilename in glob.glob(os.path.join(args.dir, args.name)):
            lvl = 1
            with open(gwasfilename, 'r') as gwasfile:
                sys.stdout.write("reading file '%s'\n" % gwasfilename)
                gwasvalues = []
                gwaspvalues = []
                linenr = 0
                for line in gwasfile:
                    linenr += 1
                    cols = re.split(',|\t', line)
                    if len(cols) == 6:
                        try:
                            gwaspvalues.append(float(cols[2]))
                            gwasvalues.append(cols)
                        except ValueError:
                            # print("could not read p-value in file: {} (line {})", gwasfilename, linenr)
                            continue
                            #                     try:
                            #                         cchr=int(cols[0])
                            #                         cpos=int(cols[1])
                            #                         csco=float(cols[2])
                            #                         cmaf=float(cols[3])
                            #                         cmac=int(cols[4])
                            #                         cvar=float(cols[5])
                            #                         gwasvalues.append([cchr,cpos,csco,cmaf,cmac,cvar])
                            #                     except ValueError:
                            #                         #sys.stderr.write("error in line %d: '%s'\n" % (linenr,line.strip()))
                            #                         continue

                            #            gwasvalues=np.array(gwasvalues)
                            #             print(gwasfilename)
                            #             print(gwasvalues[:,2])

            bonf_thres = args.fdr / float(len(gwaspvalues))
            bh_thres = mt.get_bh_thres(gwaspvalues, args.fdr)['thes_pval']
            bhy_thres = mt.get_bhy_thres(gwaspvalues, args.fdr)['thes_pval']
            try:
                used_threshold = float(args.pvalue_threshold)
                sys.stdout.write("using custom threshold: {:e}\n".format(used_threshold))
            except ValueError:
                if args.pvalue_threshold.lower() == 'bonf':
                    used_threshold = bonf_thres
                    sys.stdout.write("using bonferroni threshold (fdr {}): {:e}\n".format(args.fdr, used_threshold))
                elif args.pvalue_threshold.lower() == 'bh':
                    used_threshold = bh_thres
                    sys.stdout.write(
                        "using benjamini-hochberg threshold (fdr {}): {:e}\n".format(args.fdr, used_threshold))
                elif args.pvalue_threshold.lower() == 'bhy':
                    used_threshold = bhy_thres
                    sys.stdout.write("using benjamini-hochberg-yekutieli threshold (fdr {}): {:e}\n"
                                     .format(args.fdr,used_threshold))
                else:
                    raise Exception("unkown p-value threshold method.")

            # fdr_rejected, fdr_adjusted = sm.fdrcorrection0(gwaspvalues, alpha=0.05, method='negcorr', is_sorted=False)
            # print(len(sm.multipletests(gwaspvalues, alpha=0.05, method='fdr_by', is_sorted=False, returnsorted=False)))
            fdr_rejected, fdr_adjusted, dummy1, dummy2 = sm.multipletests(gwaspvalues, alpha=args.fdr, method='fdr_bh',
                                                                          is_sorted=False, returnsorted=False)
            sys.stdout.write("checking pvalues: \n")

            passed_cnt = 0
            for idx in range(len(gwasvalues)):
                # Lotus pvalue files now use chromosome numbering from 0:6
                # conversion is not needed anymore

                # if dbid == 'Lj'.encode():
                #     gw_chr = str(int(gwasvalues[idx][0])-1)
                # else:
                #     gw_chr = gwasvalues[idx][0]
                gw_chr = gwasvalues[idx][0]

                gw_pos = int(gwasvalues[idx][1])
                gw_pval = float(gwasvalues[idx][2])
                gw_mac = int(gwasvalues[idx][4])

                if gw_pval <= used_threshold and gw_mac >= args.minor_allele_count:
                    passed_cnt += 1
                    dbextract.flush()
                    genes = dbextract.extract_loc_uddist(gw_chr, gw_pos, args.udistance, args.ddistance)
                    # genes = dbextract.get_genes()
                    sys.stdout.write("peak: chr{}, pos {} -> {} genes in range\n".format(gw_chr, gw_pos, len(genes)))

                    if len(genes) == 0:
                        ostream.write("%s\t" % os.path.basename(gwasfilename))
                        ostream.write("%s\t" % gw_chr)
                        ostream.write("%s\t" % gw_pos)
                        ostream.write("%s\t" % gw_pval)
                        ostream.write("%s\t" % fdr_rejected[idx])
                        ostream.write("%f\t" % fdr_adjusted[idx])
                        ostream.write("%e\t" % bonf_thres)
                        ostream.write("%e\t" % bh_thres)
                        ostream.write("%e\t" % bhy_thres)
                        ostream.write("NA\t")
                        ostream.write("NA\t")
                        ostream.write("?\t")
                        ostream.write("NA\t")
                        ostream.write("NoGeneFound\t")
                        ostream.write("NoGeneFound\t")
                        ostream.write("NoGeneFound\t")
                        # ostream.write("%s\t" % gene.shortsym)
                        # ostream.write("%s\t" % gene.longsym)
                        ostream.write("None\t")
                        ostream.write("None\t")
                        ostream.write("NoGeneFound\n")
                        ostream.flush()
                        continue
                    for gene in genes:
                        ostream.write("%s\t" % os.path.basename(gwasfilename))
                        ostream.write("%s\t" % gw_chr)
                        ostream.write("%s\t" % gw_pos)
                        ostream.write("%s\t" % gw_pval)
                        ostream.write("%s\t" % fdr_rejected[idx])
                        ostream.write("%f\t" % fdr_adjusted[idx])
                        ostream.write("%e\t" % bonf_thres)
                        ostream.write("%e\t" % bh_thres)
                        ostream.write("%e\t" % bhy_thres)
                        ostream.write("%d\t" % gene.start)
                        ostream.write("%d\t" % gene.end)
                        ostream.write("%s\t" % gene.strand)
                        if gene.strand == '+':
                            ostream.write("%d\t" % (gene.start - gw_pos))
                        else:
                            ostream.write("%d\t" % (gw_pos - gene.start))

                        if gene.start <= gw_pos <= gene.end:
                            relpos = "in gene"
                        elif gw_pos < gene.start:
                            if gene.strand == '+':
                                relpos = "upstream"
                            else:
                                relpos = "downstream"
                        else:
                            if gene.strand == '+':
                                relpos = "downstream"
                            else:
                                relpos = "upstream"
                        ostream.write("%s\t" % relpos)
                        ostream.write("%s\t" % gene.id)
                        ostream.write("%s\t" % gene.feature)
                        # ostream.write("%s\t" % gene.shortsym)
                        # ostream.write("%s\t" % gene.longsym)
                        ostream.write("%s\t" % gene.sequencetype)
                        ostream.write("None\t")
                        ostream.write("%s\n" % gene.attribute)

                        if args.depth >= 1:
                            for rna in gene.rna:
                                ostream.write("%s\t" % os.path.basename(gwasfilename))
                                ostream.write("%s\t" % gw_chr)
                                ostream.write("%s\t" % gw_pos)
                                ostream.write("%s\t" % gw_pval)
                                ostream.write("%s\t" % fdr_rejected[idx])
                                ostream.write("%f\t" % fdr_adjusted[idx])
                                ostream.write("%e\t" % bonf_thres)
                                ostream.write("%e\t" % bh_thres)
                                ostream.write("%e\t" % bhy_thres)
                                ostream.write("%d\t" % rna.start)
                                ostream.write("%d\t" % rna.end)
                                ostream.write("%s\t" % rna.strand)
                                if rna.strand == '+':
                                    ostream.write("%d\t" % (rna.start - gw_pos))
                                else:
                                    ostream.write("%d\t" % (gw_pos - rna.end))

                                if rna.start <= gw_pos <= rna.end:
                                    relpos = "in gene"
                                elif gw_pos < rna.start:
                                    if rna.strand == '+':
                                        relpos = "upstream"
                                    else:
                                        relpos = "downstream"
                                else:
                                    if rna.strand == '+':
                                        relpos = "downstream"
                                    else:
                                        relpos = "upstream"
                                ostream.write("%s\t" % relpos)
                                ostream.write("%s\t" % rna.id)
                                ostream.write("%s\t" % rna.feature)
                                # ostream.write("%s\t" % gene.shortsym)
                                ostream.write("%s\t" % rna.sequencetype)
                                ostream.write("%s\t" % rna.short_annotation)
                                ostream.write("%s\n" % rna.attribute)
                                ostream.flush()
            sys.stdout.write("{} peaks passed.\n".format(passed_cnt))
        if isfile:
            ostream.close()

    @staticmethod
    def remap_hdf5_results(args):
        for h5filepath in glob.glob(os.path.join(args.dir, args.name)):
            GwasData.remap_hdf5_results(h5filepath, args.mapfile)

    def parse_arguments(self):
        mainparser = argparse.ArgumentParser(description='tair database suite',
                                             prog="genehunter")  # ,usage='genehunter command [<args>]')
        subparsers = mainparser.add_subparsers(dest='command', help='subcommand help')
        subparsers.required = True

        # createdb
        createparser = subparsers.add_parser('createdb', help='create new database from tair10 files')
        createparser.add_argument('--gff', required=True, help='path to gff file')
        createparser.add_argument('--desc', help='path to functional descriptions file')
        createparser.add_argument('--aliases', help='path to gene aliases file')
        createparser.add_argument('--sorf', help='path to sORF file')
        createparser.add_argument('-a', '--add', action="store_true", help='add entries to existing database')
        createparser.add_argument('-o', '--output', required=True, help='path to output file')
        createparser.set_defaults(func=self.create_db)

        extract_by_loc = subparsers.add_parser('extractloc', help='extract elements by locus')
        extract_by_loc.add_argument('--db', required=True, help='path to database')
        extract_by_loc.add_argument('--depth', type=int, default=1,
                                    help='defines what to print. 0=genes only, 1=genes and rna, 2=all features')
        extract_by_loc.add_argument('--loc1', type=int,
                                    help='locus 1 (has different meanings depending on other options)')
        extract_by_loc.add_argument('--loc2', type=int, help='locus 2 (different meanings depending on other options)')
        extract_by_loc.add_argument('--chr', help='chromosome')
        extract_by_loc.add_argument('-c', action="store_true", help='centered around loc1, interval = +/-loc2')
        extract_by_loc.add_argument('-i', action="store_true", help='from locus interval. loc1=start, loc2=end')
        extract_by_loc.add_argument('--file',
                                    help='read positions from file. The file should contain lines with 4 columns each. \
                                    (chromosome, position, upstream interval,downstream interval). If only 2 columns \
                                    are present, the interval will be pos-loc1 to pos+loc2')
        extract_by_loc.add_argument('-o', '--output', help='path to output file. Will print to stdout if omitted')
        extract_by_loc.set_defaults(func=self.extract_loc)

        extract_by_agi = subparsers.add_parser('extractagi', help='extract elements by AGI')
        extract_by_agi.add_argument('--db', required=True, help='path to database')
        extract_by_agi.add_argument('--depth', type=int, default=1,
                                    help='defines what to print. 0=genes only, 1=genes and rna, 2=all features')
        extract_by_agi.add_argument('--agi', help='single AGI')
        extract_by_agi.add_argument('--file', help='read AGIs from file. The file must contain one AGI per line.')
        extract_by_agi.add_argument('-o', '--output', help='path to output file. Will print to stdout if omitted')
        extract_by_agi.set_defaults(func=self.extract_agi)

        hunterparser = subparsers.add_parser('hunt', help='run gene hunter on pvals and extract gene information')
        hunterparser.add_argument('--db', required=True, help='path to tair10 database')
        hunterparser.add_argument('--dir', required=True, help='directory where to look for pval files')
        hunterparser.add_argument('--name', default="*.pvals",
                                  help='name identifier for the pval files. Unix like globs are allow (e.g. "name*")')
        hunterparser.add_argument('--depth', type=int, default=1,
                                  help='defines what to print. 0=genes only, 1=genes and rna, 2=all features')
        hunterparser.add_argument('-u', '--udistance', type=int, default=4000,
                                  help='maximal upstream distance from TSS (default=4000)')
        hunterparser.add_argument('-d', '--ddistance', type=int, default=4000,
                                  help='maximal downstream distance from TSS (default=4000)')
        hunterparser.add_argument('-P', '--pvalue_threshold', default='1.0e-6', type=float,
                                  help='SNP p-value threshold (default=1.0e-6).') # If the given argument is \"bonf\", \
                                  # \"bh\" or \"bhy\", only values below the calculated threshold are included.')
        hunterparser.add_argument('-M', '--minor_allele_count', type=int, default=10,
                                  help='minor allele count threshold (default=10)')
        hunterparser.add_argument('-o', '--output', default=None, help='Prefix for output file. Will print to stdout if omitted.')
        hunterparser.add_argument('-f', '--fdr', type=float, default=0.05,
                                  help="alpha for fdr threshold calculation")
        hunterparser.set_defaults(func=self.extract_hunter)

        remapparser = subparsers.add_parser('remap_chrs', help='remap integer chromosome names to real names')
        remapparser.add_argument('-d', '--dir', required=True, help='directory where to look for pval files')
        remapparser.add_argument('-n', '--name', default="*.hdf5", required=True, help='hdf5 file(s) to be changed')
        remapparser.add_argument('-m', '--mapfile', required=True, help='comma separated txt file containing the mapping')
        remapparser.set_defaults(func=self.remap_hdf5_results)

        statusparser = subparsers.add_parser('stats', help='get some statistics about database')
        statusparser.add_argument('--db', required=True, help='path to gene annotation database')
        statusparser.set_defaults(func=self.print_stats)

        args = mainparser.parse_args()
        args.func(args)


# if __name__ == '__main__':
#     extractor = GeneAnnotationDbExtractor("/home/GMI/christian.goeschl/devel/pycharm/GeneHunter/db/Mt4.0v2_20140818_1100.sqlite")
#     extractor.extract_test()

if __name__ == '__main__':
    main()
