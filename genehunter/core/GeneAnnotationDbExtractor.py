# -*- coding: utf-8 -*-

"""genehunter.tairdb_extractor: module for querying the database."""

import os
import sys
import re
import unicodedata
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.exc import *
# from sqlalchemy import and_, or_
from genehunter.core.HunterData import OutputData
from .GeneAnnotationDbModels import Gene


def _re_fn(expr, item):
    reg = re.compile(expr, re.I)
    return reg.search(item) is not None


class GeneAnnotationDbExtractor:
    def __init__(self, dbpath):
        self.dbpath = dbpath
        self.genes = []
        self.engine = create_engine('sqlite:///' + dbpath, echo=False)
        if not os.path.isfile(dbpath):
            sys.stderr.write("could not open database %s" % dbpath)
            return

        conn = self.engine.connect()
        conn.connection.create_function('regexp', 2, _re_fn)
        session = sessionmaker()
        session.configure(bind=conn)  # (bind=self.engine)
        self.session = session()

    def flush(self):
        self.genes = []

    def get_genes(self):
        return self.genes

    def print_stats(self):
        allgenes = self.session.query(Gene.id).all()
        sys.stdout.write("database: {} contains:\n".format(self.dbpath))
        sys.stdout.write("\t {:d} genes in total.\n".format(len(allgenes)))

        ljcnt = 0
        ljchlorocnt = 0
        ljmitocnt = 0
        othercnt = 0
        ljPattern = '^Lj\dg'
        ljchloroPattern = '^Ljchlorog'
        ljmitoPattern = '^Ljmitog'
        ljset = []
        for id in allgenes:
            if re.match(ljPattern,id[0]):
                ljcnt += 1
                ljset.append(id[0])
                if len(id[0]) != 13:
                    pass
            elif re.match(ljchloroPattern,id[0]):
                ljchlorocnt += 1
            elif re.match(ljmitoPattern, id[0]):
                ljmitocnt += 1
            else:
                othercnt += 1
        sys.stdout.write("\t {:d} Lj... genes.\n".format(ljcnt))
        sys.stdout.write("\t {:d} unique Lj... genes.\n".format(len(set(ljset))))
        sys.stdout.write("\t {:d} Ljchloro... genes.\n".format(ljchlorocnt))
        sys.stdout.write("\t {:d} Ljmito... genes.\n".format(ljmitocnt))
        sys.stdout.write("\t {:d} other genes.\n".format(othercnt))

    def get_database_id(self):
        geneid = self.session.query(Gene.id).distinct().first()[0]
        dbid = unicodedata.normalize('NFKD', geneid).encode('ascii', 'ignore')[0:2]
        return dbid

    def extract_by_agi(self, agi):
        eagi = agi.split('.')[0]
        #         if(len(eagi)>1):
        #             try:
        #                 emodel=int(eagi[1])
        #             except ValueError:
        #                 emodel=None
        self.genes.extend(self.session.query(Gene).filter(Gene.id == eagi).all())

    def extract_by_agi_list(self, agilist):
        for item in agilist:
            self.extract_by_agi(item)

    def extract_by_loc(self, eloc_start, eloc_end, echrom):
        try:
            int(echrom)
            chr_int_search = True
        except ValueError:
            chr_int_search = False

        if chr_int_search:
            regexstr = "[a-z_]{0,}" + echrom
            self.genes.extend(self.session.query(Gene).filter(Gene.seqname.op('regexp')(regexstr)).filter(
                Gene.start.between(eloc_start, eloc_end) |
                Gene.end.between(eloc_start, eloc_end)))
        else:
            self.genes.extend(self.session.query(Gene).filter(Gene.seqname == echrom).filter(
                Gene.start.between(eloc_start, eloc_end) |
                Gene.end.between(eloc_start, eloc_end)))

            # ((Gene.end >= eloc_start) & (Gene.end <= eloc_end)) |
            # ((Gene.start >= eloc_start) & (Gene.start <= eloc_end)) |
            # ((Gene.start < eloc_start) & (Gene.end > eloc_end))))

            #     self.genes.extend(self.session.query(TairGene).filter(TairGene.chromosome.op('regexp')('.*3'))
            #                       .filter(((TairGene.loc_end >= eloc_start) & (TairGene.loc_end <= eloc_end)) |
            #                               ((TairGene.loc_start >= eloc_start) & (TairGene.loc_start <= eloc_end)) |
            #                               ((TairGene.loc_start < eloc_start) & (TairGene.loc_end > eloc_end))))
            #
            # else:
            #     self.genes.extend(self.session.query(TairGene).filter(TairGene.chromosome == echrom)
            #                       .filter(((TairGene.loc_end >= eloc_start) & (TairGene.loc_end <= eloc_end)) |
            #                               ((TairGene.loc_start >= eloc_start) & (TairGene.loc_start <= eloc_end)) |
            #                               ((TairGene.loc_start < eloc_start) & (TairGene.loc_end > eloc_end))))

    # def extract_oriented_loc(self, eloc_start, eloc_end, echrom, eorient):
    #     self.genes.extend(self.session.query(TairGene).filter(TairGene.chromosome == echrom)
    #                       .filter(TairGene.orientation == eorient)
    #                       .filter(((TairGene.loc_end >= eloc_start) & (TairGene.loc_end <= eloc_end)) |
    #                               ((TairGene.loc_start >= eloc_start) & (TairGene.loc_start <= eloc_end)) |
    #                               ((TairGene.loc_start < eloc_start) & (TairGene.loc_end > eloc_end))))
    #
    def extract_by_loc(self, inputdata, infodepth=1):
        output_df = OutputData()
        for ix in inputdata.df.index:
            dseries = inputdata.df.loc[ix]
            genes = self.extract_loc_uddist(dseries["Chromosome"], dseries["SNP_pos"], dseries["uDist"], dseries["dDist"])
            for gene in genes:
                out_series = pd.Series()
                out_series.iloc[0:3] = dseries.iloc[0:3]
                out_series["Gene_start"] = gene.start
                out_series["Gene_end"] = gene.end
                out_series["Gene_orientation"] = gene.strand
                if gene.strand == '+':
                    out_series["Relative_Distance"] = dseries["SNP_pos"] - gene.start
                else:
                    out_series["Relative_Distance"] = gene.start - dseries["SNP_pos"]
                if gene.start <= dseries["SNP_pos"] <= gene.end:
                    relpos = "in gene"
                elif dseries["SNP_pos"] < gene.start:
                    if gene.strand == '+':
                        relpos = "upstream"
                    else:
                        relpos = "downstream"
                else:
                    if gene.strand == '+':
                        relpos = "downstream"
                    else:
                        relpos = "upstream"
                out_series["SNP_relative_position"] = relpos
                out_series["target_AGI"] = gene.id
                out_series["target_element_type"] = gene.feature
                out_series["target_sequence_type"] = gene.sequencetype
                out_series["target_attributes"] = gene.attribute
                output_df.add_dataset(out_series)

                if infodepth > 0:
                    rna = gene.rna
                    out_series["Gene_start"] = rna.start
                    out_series["Gene_end"] = rna.end
                    out_series["Gene_orientation"] = rna.strand
                    if rna.strand == '+':
                        out_series["Relative_Distance"] = dseries["SNP_pos"] - rna.start
                    else:
                        out_series["Relative_Distance"] = rna.start - dseries["SNP_pos"]
                    if rna.start <= dseries["SNP_pos"] <= rna.end:
                        relpos = "in gene"
                    elif dseries["SNP_pos"] < rna.start:
                        if rna.strand == '+':
                            relpos = "upstream"
                        else:
                            relpos = "downstream"
                    else:
                        if rna.strand == '+':
                            relpos = "downstream"
                        else:
                            relpos = "upstream"
                    out_series["SNP_relative_position"] = relpos
                    out_series["target_AGI"] = rna.id
                    out_series["target_element_type"] = rna.feature
                    out_series["target_sequence_type"] = rna.sequencetype
                    out_series["target_attributes"] = rna.attribute
                    output_df.add_dataset(out_series)
                # ostream.write(str(gene.seqname) + "\t")
                # ostream.write(str(gene.start) + "\t")
                # ostream.write(str(gene.end) + "\t")
                # ostream.write(str(gene.strand) + "\t")
                # ostream.write(str(gene.id) + "\t")
                # ostream.write(str(gene.feature) + "\t")
                # # ostream.write(str(gene.shortsym) + "\t")
                # # ostream.write(str(gene.longsym) + "\t")
                # ostream.write("\t")
                # ostream.write(str(gene.attribute) + "\n")
        return output_df

    def extract_loc_uddist(self, echr, eloc, udist, ddist):
        genes = []
        try:
            int(echr)
            chr_int_search = True
        except ValueError:
            chr_int_search = False

        if chr_int_search:
            regexstr = "[a-z_]{0,}" + echr
            genes.extend(self.session.query(Gene).filter(Gene.seqname.op('regexp')(regexstr))
                         .filter(Gene.strand == '+')
                         .filter(Gene.start.between(eloc - udist, eloc + ddist) |
                                 Gene.end.between(eloc - udist, eloc + ddist) |
                                 ((Gene.start <= eloc) & (Gene.end >= eloc))))
            genes.extend(self.session.query(Gene).filter(Gene.seqname.op('regexp')(regexstr))
                              .filter(Gene.strand == '-')
                              .filter(Gene.start.between(eloc - ddist, eloc + udist) |
                                      Gene.end.between(eloc - ddist, eloc + udist) |
                                      ((Gene.start <= eloc) & (Gene.end >= eloc))))
        else:
            genes.extend(self.session.query(Gene).filter(Gene.seqname == echr)
                              .filter(Gene.strand == '+')
                              .filter(Gene.start.between(eloc-udist, eloc+ddist) |
                                      Gene.end.between(eloc-udist, eloc+ddist) |
                                      ((Gene.start <= eloc) & (Gene.end >= eloc))))
            genes.extend(self.session.query(Gene).filter(Gene.seqname == echr)
                              .filter(Gene.strand == '-')
                              .filter(Gene.start.between(eloc - ddist, eloc + udist) |
                                      Gene.end.between(eloc - ddist, eloc + udist) |
                                      ((Gene.start <= eloc) & (Gene.end >= eloc))))

        # self.genes.extend(self.session.query(TairGene)
        #                   .filter(TairGene.chromosome == echr)
        #                   .filter(TairGene.orientation == '+')
        #                   .filter(
        #     ((TairGene.loc_start - udist) <= eloc) & ((TairGene.loc_end + ddist) >= eloc)).all())
        # self.genes.extend(self.session.query(TairGene)
        #                   .filter(TairGene.chromosome == echr)
        #                   .filter(TairGene.orientation == '-')
        #                   .filter(
        #     ((TairGene.loc_end + udist) >= eloc) & ((TairGene.loc_start - ddist) <= eloc)).all())
        return genes

    def write_results(self, outpath=None, depth=0):
        if outpath is not None:
            ostream = open(outpath, 'w')
            isfile = True
        else:
            ostream = sys.stdout
            isfile = False

        ostream.write("Chromosome\tSNP_pos\tLoc_start\tLoc_end\tOrientation\tAGI\tType\tAttributes\n")
        # ostream.write("Short_sym\tLong_sym\tShort_Desc\tLong_Desc\n")
        for gene in self.genes:
            ostream.write(str(gene.seqname) + "\t")
            ostream.write(str(gene.start) + "\t")
            ostream.write(str(gene.end) + "\t")
            ostream.write(str(gene.strand) + "\t")
            ostream.write(str(gene.id) + "\t")
            ostream.write(str(gene.feature) + "\t")
            # ostream.write(str(gene.shortsym) + "\t")
            # ostream.write(str(gene.longsym) + "\t")
            ostream.write("\t")
            ostream.write(str(gene.attribute) + "\n")

            if depth >= 1:
                for rna in gene.rna:
                    # desc = rna.description
                    ostream.write(str(rna.seqname) + "\t")
                    ostream.write(str(rna.start) + "\t")
                    ostream.write(str(rna.end) + "\t")
                    ostream.write(str(rna.strand) + "\t")
                    ostream.write(str(rna.id) + "\t")
                    ostream.write(str(rna.feature) + "\t")
                    # ostream.write(str(rna.shortsym) + "\t")
                    # ostream.write(str(rna.longsym) + "\t")
                    # ostream.write(str(desc.shortdesc) + "\t")
                    ostream.write(str(rna.attribute) + "\n")

                    if depth >= 2:
                        for feature in rna.features:
                            ostream.write(str(feature.seqname) + "\t")
                            ostream.write(str(feature.start) + "\t")
                            ostream.write(str(feature.end) + "\t")
                            ostream.write(str(feature.strand) + "\t")
                            ostream.write(str(feature.id) + "\t")
                            ostream.write(str(feature.feature) + "\t")
                            # ostream.write(str(gene.shortsym) + "\t")
                            # ostream.write(str(gene.longsym) + "\t")
                            # ostream.write(str(desc.shortdesc) + "\t")
                            ostream.write(str(feature.attribute) + "\n")
        if isfile:
            ostream.close()
