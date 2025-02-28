# -*- coding: utf-8 -*-

"""tairdbsuite.tairdbsuite: module doing the argument parsing and dispatching the commands."""

import sys
import os
import re
import glob
import argparse
import statsmodels.sandbox.stats.multicomp as sm
# import numpy as np
import tairdbsuite.core.mtcorr as mt
from tairdbsuite.core.TairDBCreator import TairDBCreator
from tairdbsuite.core.TairDBExtractor import TairDBExtractor

__version__ = "0.1.2"


def main():
    sys.stdout.write("tair database suite - version %s\n" % __version__)
    TairdbSuite()


class TairdbSuite(object):
    def __init__(self):
        self.parse_arguments()

    @staticmethod
    def create_db(args):
        if os.path.exists(args.output):
            raise Exception("database file '%s' exists!\nChoose another output file. Terminating!\n" % args.output)
        if args.gff is not None and not os.path.exists(args.gff):
            raise Exception("cannot open GFF file '%s'! Terminating!\n" % args.gff)
        if args.desc is not None and not os.path.exists(args.desc):
            raise Exception("cannot open descriptions file '%s'! Terminating!\n" % args.desc)
        if args.aliases is not None and not os.path.exists(args.aliases):
            raise Exception("cannot open aliases file '%s'! Terminating!\n" % args.aliases)
        if args.sorf is not None and not os.path.exists(args.sorf):
            raise Exception("cannot open sORF file '%s'! Terminating!\n" % args.sorf)
        creator = TairDBCreator(args.output)
        creator.create_gff_entries(args.gff)
        if args.desc is not None:
            creator.create_desc_entries(args.desc)
        if args.aliases is not None:
            creator.create_name_entries(args.aliases)
        if args.sorf is not None:
            creator.create_sorf_entries(args.sorf)
        return

    @staticmethod
    def print_stats(args):
        if not os.path.exists(args.db):
            raise Exception("cannot open database '%s'!" % args.db)

        extractor = TairDBExtractor(args.db)
        extractor.print_stats()

    @staticmethod
    def extract_loc(args):
        intervals = []
        if not os.path.exists(args.db):
            raise Exception("cannot open database '%s'!" % args.db)

        if args.file is not None:
            with open(args.file, 'r') as ifile:
                for line in ifile:
                    cols = re.split(',|\t', line)
                    chrom = cols[0].strip()
                    pos = int(cols[1].strip())
                    if len(cols) == 2:
                        intervals.append((chrom, pos - args.loc1, pos + args.loc2))
                    else:
                        intervals.append((chrom, pos - int(cols[2].strip()), pos + int(cols[3].strip())))
        elif args.c:
            intervals.append((args.chr, args.loc1 - args.loc2, args.loc1 + args.loc2))
        elif args.i:
            intervals.append((args.chr, args.loc1, args.loc2))

        extractor = TairDBExtractor(args.db)
        for interval in intervals:
            extractor.extract_by_loc(interval[1], interval[2], interval[0])
        extractor.write_results(args.output, 1)

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

        extractor = TairDBExtractor(args.db)
        for agi in agis:
            extractor.extract_by_agi(agi)
        extractor.write_results(args.output, 1)
        return

    @staticmethod
    def extract_hunter(args):
        dbextract = TairDBExtractor(args.db)

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
        ostream.write("short_symbol\tlong_symbol\tshort_description\tlong_description\n")

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

            bonf_thres = args.fdr/float(len(gwaspvalues))
            bh_thres = mt.get_bh_thres(gwaspvalues, args.fdr)['thes_pval']
            bhy_thres = mt.get_bhy_thres(gwaspvalues, args.fdr)['thes_pval']
            try:
                used_threshold = float(args.pvalue_threshold)
                sys.stdout.write("using custom threshold: {:e}\n".format(used_threshold))
            except ValueError:
                if args.pvalue_threshold.lower() == 'bonf':
                    used_threshold = bonf_thres
                    sys.stdout.write("using bonferroni threshold (fdr {}): {:e}\n".format(args.fdr,used_threshold))
                elif args.pvalue_threshold.lower() == 'bh':
                    used_threshold = bh_thres
                    sys.stdout.write("using benjamini-hochberg threshold (fdr {}): {:e}\n".format(args.fdr, used_threshold))
                elif args.pvalue_threshold.lower() == 'bhy':
                    used_threshold = bhy_thres
                    sys.stdout.write("using benjamini-hochberg-yekutieli threshold (fdr {}): {:e}\n".format(args.fdr, used_threshold))
                else:
                    raise Exception("unkown p-value threshold method.")

            # fdr_rejected, fdr_adjusted = sm.fdrcorrection0(gwaspvalues, alpha=0.05, method='negcorr', is_sorted=False)
            # print(len(sm.multipletests(gwaspvalues, alpha=0.05, method='fdr_by', is_sorted=False, returnsorted=False)))
            fdr_rejected, fdr_adjusted, dummy1, dummy2 = sm.multipletests(gwaspvalues, alpha=args.fdr, method='fdr_bh',
                                                                          is_sorted=False, returnsorted=False)
            sys.stdout.write("checking pvalues: \n")

            passed_cnt = 0
            for idx in range(len(gwasvalues)):
                gw_chr = gwasvalues[idx][0]
                gw_pos = int(gwasvalues[idx][1])
                gw_pval = float(gwasvalues[idx][2])
                gw_mac = int(gwasvalues[idx][4])

                if gw_pval <= used_threshold and gw_mac >= args.minor_allele_count:
                    passed_cnt += 1
                    dbextract.flush()
                    dbextract.extract_loc_uddist(gw_chr, gw_pos, args.udistance, args.ddistance)
                    genes = dbextract.get_genes()
                    sys.stdout.write("peak: chr{}, pos {} -> {} genes in range\n".format(gw_chr, gw_pos, len(genes)))

                    for gene in genes:
                        ostream.write("%s\t" % os.path.basename(gwasfilename))
                        ostream.write("%s\t" % gw_chr)
                        ostream.write("%s\t" % gw_pos)
                        ostream.write("%s\t" % gw_pval)
                        ostream.write("%f\t" % fdr_rejected[idx])
                        ostream.write("%f\t" % fdr_adjusted[idx])
                        ostream.write("%e\t" % bonf_thres)
                        ostream.write("%e\t" % bh_thres)
                        ostream.write("%e\t" % bhy_thres)
                        ostream.write("%d\t" % gene.loc_start)
                        ostream.write("%d\t" % gene.loc_end)
                        ostream.write("%s\t" % gene.orientation)
                        if gene.orientation == '+':
                            ostream.write("%d\t" % abs(gene.loc_start - gw_pos))
                        else:
                            ostream.write("%d\t" % abs(gene.loc_end - gw_pos))

                        if gene.loc_start <= gw_pos <= gene.loc_end:
                            relpos = "in gene"
                        elif gw_pos < gene.loc_start:
                            if gene.orientation == '+':
                                relpos = "upstream"
                            else:
                                relpos = "downstream"
                        else:
                            if gene.orientation == '+':
                                relpos = "downstream"
                            else:
                                relpos = "upstream"
                        ostream.write("%s\t" % relpos)
                        ostream.write("%s\t" % gene.agi)
                        ostream.write("%s\t" % gene.type)
                        ostream.write("%s\t" % gene.shortsym)
                        ostream.write("%s\t" % gene.longsym)
                        ostream.write("None\t")
                        ostream.write("%s\n" % gene.attributes)

                        if lvl >= 1:
                            for child in gene.children:
                                ostream.write("%s\t" % os.path.basename(gwasfilename))
                                ostream.write("%s\t" % gw_chr)
                                ostream.write("%s\t" % gw_pos)
                                ostream.write("%s\t" % gw_pval)
                                ostream.write("%f\t" % fdr_rejected[idx])
                                ostream.write("%f\t" % fdr_adjusted[idx])
                                ostream.write("%e\t" % bonf_thres)
                                ostream.write("%e\t" % bh_thres)
                                ostream.write("%e\t" % bhy_thres)
                                ostream.write("%d\t" % child.loc_start)
                                ostream.write("%d\t" % child.loc_end)
                                ostream.write("%s\t" % gene.orientation)
                                if gene.orientation == '+':
                                    ostream.write("%d\t" % abs(child.loc_start - gw_pos))
                                else:
                                    ostream.write("%d\t" % abs(child.loc_end - gw_pos))

                                if child.loc_start <= gw_pos <= child.loc_end:
                                    relpos = "in gene"
                                elif gw_pos < gene.loc_start:
                                    if gene.orientation == '+':
                                        relpos = "upstream"
                                    else:
                                        relpos = "downstream"
                                else:
                                    if gene.orientation == '+':
                                        relpos = "downstream"
                                    else:
                                        relpos = "upstream"
                                ostream.write("%s\t" % relpos)
                                ostream.write("%s.%d\t" % (gene.agi, child.model))
                                ostream.write("%s\t" % child.type)
                                ostream.write("%s\t" % gene.shortsym)
                                ostream.write("%s\t" % gene.longsym)
                                ostream.write("%s\t" % child.description.shortdesc)
                                ostream.write("%s\n" % child.description.longdesc)
                                ostream.flush()
            sys.stdout.write("{} peaks passed.\n".format(passed_cnt))
        if isfile:
            ostream.close()

    def parse_arguments(self):
        mainparser = argparse.ArgumentParser(description='tair database suite',
                                             prog="tairdbsuite")  # ,usage='tairdbsuite command [<args>]')
        subparsers = mainparser.add_subparsers(dest='command', help='subcommand help')
        subparsers.required = True

        # createdb
        createparser = subparsers.add_parser('createdb', help='create new database from tair10 files')
        createparser.add_argument('--gff', required=True, help='path to gff file')
        createparser.add_argument('--desc', help='path to functional descriptions file')
        createparser.add_argument('--aliases', help='path to gene aliases file')
        createparser.add_argument('--sorf', help='path to sORF file')
        createparser.add_argument('-o', '--output', required=True, help='path to output file')
        createparser.set_defaults(func=self.create_db)

        extract_by_loc = subparsers.add_parser('extractloc', help='extract elements by locus')
        extract_by_loc.add_argument('--db', required=True, help='path to database')
        extract_by_loc.add_argument('--loc1', type=int,
                                    help='locus 1 (has different meanings depending on other options)')
        extract_by_loc.add_argument('--loc2', type=int, help='locus 2 (different meanings depending on other options)')
        extract_by_loc.add_argument('--chr', type=int, help='chromosome')
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
        extract_by_agi.add_argument('--agi', help='single AGI')
        extract_by_agi.add_argument('--file', help='read AGIs from file. The file must contain one AGI per line.')
        extract_by_agi.add_argument('-o', '--output', help='path to output file. Will print to stdout if omitted')
        extract_by_agi.set_defaults(func=self.extract_agi)

        hunterparser = subparsers.add_parser('hunter', help='run gene hunter on pvals and extract gene information')
        hunterparser.add_argument('--db', required=True, help='path to tair10 database')
        hunterparser.add_argument('--dir', required=True, help='directory where to look for pval files')
        hunterparser.add_argument('--name', default="*.pvals",
                                  help='name identifier for the pval files. Unix like globs are allow (e.g. "name*")')
        hunterparser.add_argument('-u', '--udistance', type=int, default=4000,
                                  help='maximal upstream distance from TSS (default=4000)')
        hunterparser.add_argument('-d', '--ddistance', type=int, default=4000,
                                  help='maximal downstream distance from TSS (default=4000)')
        hunterparser.add_argument('-P', '--pvalue_threshold', default='1.0e-6',
                                  help='SNP p-value threshold (default=1.0e-6). If the given argument is \"bonf\", \
                                  \"bh\" or \"bhy\", only values below the calculated threshold are included.')
        hunterparser.add_argument('-M', '--minor_allele_count', type=int, default=10,
                                  help='minor allele count threshold (default=10)')
        hunterparser.add_argument('-o', '--output', help='path to output file. Will print to stdout if omitted.')
        hunterparser.add_argument('-f', '--fdr', type=float, default=0.05,
                                  help="alpha for fdr threshold calculation")
        hunterparser.set_defaults(func=self.extract_hunter)

        statusparser = subparsers.add_parser('stats', help='get some statistics about database')
        statusparser.add_argument('--db', required=True, help='path to gene annotation database')
        statusparser.set_defaults(func=self.print_stats)

        args = mainparser.parse_args()
        args.func(args)


if __name__ == '__main__':
    main()
