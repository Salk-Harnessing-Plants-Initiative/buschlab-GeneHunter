# -*- coding: utf-8 -*-

"""tairdbsuite.tairdb_extractor: module for querying the database."""

import os
import sys
import re
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.exc import *
# from sqlalchemy import and_, or_

from .TairDBModels import TairGene  # , Tair_Level1, Tair_Level2, Tair_Desc
from .GeneAnnotationDbModels import Gene


def _re_fn(expr, item):
    reg = re.compile(expr, re.I)
    return reg.search(item) is not None


class GeneAnnotationDbExtractor:
    def __init__(self, dbpath):
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

    def extract_by_agi(self, agi):
        eagi = agi.split('.')[0]
        #         if(len(eagi)>1):
        #             try:
        #                 emodel=int(eagi[1])
        #             except ValueError:
        #                 emodel=None
        self.genes.extend(self.session.query(TairGene).filter(TairGene.agi == eagi).all())

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

    def extract_oriented_loc(self, eloc_start, eloc_end, echrom, eorient):
        self.genes.extend(self.session.query(TairGene).filter(TairGene.chromosome == echrom)
                          .filter(TairGene.orientation == eorient)
                          .filter(((TairGene.loc_end >= eloc_start) & (TairGene.loc_end <= eloc_end)) |
                                  ((TairGene.loc_start >= eloc_start) & (TairGene.loc_start <= eloc_end)) |
                                  ((TairGene.loc_start < eloc_start) & (TairGene.loc_end > eloc_end))))

    def extract_loc_uddist(self, echr, eloc, udist, ddist):
        self.genes.extend(self.session.query(TairGene)
                          .filter(TairGene.chromosome == echr)
                          .filter(TairGene.orientation == '+')
                          .filter(
            ((TairGene.loc_start - udist) <= eloc) & ((TairGene.loc_end + ddist) >= eloc)).all())
        self.genes.extend(self.session.query(TairGene)
                          .filter(TairGene.chromosome == echr)
                          .filter(TairGene.orientation == '-')
                          .filter(
            ((TairGene.loc_end + udist) >= eloc) & ((TairGene.loc_start - ddist) <= eloc)).all())

    def write_results(self, outpath=None, lvl=""):
        if outpath is not None:
            ostream = open(outpath, 'w')
            isfile = True
        else:
            ostream = sys.stdout
            isfile = False

        ostream.write("Chromosome\tLoc_start\tLoc_end\tOrientation\tAGI\tType\t")
        ostream.write("Short_sym\tLong_sym\tShort_Desc\tLong_Desc\n")
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

            # if lvl >= 1:
            #     for child in gene.children:
            #         desc = child.description
            #         ostream.write(str(gene.chromosome) + "\t")
            #         ostream.write(str(child.loc_start) + "\t")
            #         ostream.write(str(child.loc_end) + "\t")
            #         ostream.write(str(gene.orientation) + "\t")
            #         ostream.write(str(gene.agi) + "." + str(child.model) + "\t")
            #         ostream.write(str(child.type) + "\t")
            #         ostream.write(str(gene.shortsym) + "\t")
            #         ostream.write(str(gene.longsym) + "\t")
            #         ostream.write(str(desc.shortdesc) + "\t")
            #         ostream.write(str(desc.longdesc) + "\n")
            #
            #         if lvl >= 2:
            #             for subchild in child.children:
            #                 ostream.write(str(gene.chromosome) + "\t")
            #                 ostream.write(str(subchild.loc_start) + "\t")
            #                 ostream.write(str(subchild.loc_end) + "\t")
            #                 ostream.write(str(gene.orientation) + "\t")
            #                 ostream.write(str(gene.agi) + "." + str(child.model) + "\t")
            #                 ostream.write(str(subchild.type) + "\t")
            #                 ostream.write(str(gene.shortsym) + "\t")
            #                 ostream.write(str(gene.longsym) + "\t")
            #                 ostream.write(str(desc.shortdesc) + "\t")
            #                 ostream.write(str(desc.longdesc) + "\n")
        if isfile:
            ostream.close()
