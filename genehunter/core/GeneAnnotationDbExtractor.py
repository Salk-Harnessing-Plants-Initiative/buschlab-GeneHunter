# -*- coding: utf-8 -*-

"""genehunter.tairdb_extractor: module for querying the database."""

import os
import sys
import re
import unicodedata
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.exc import *
# from sqlalchemy import and_, or_

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
        sys.stdout.write("\t {:d} genes.\n".format(len(allgenes)))

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
    def extract_loc_uddist(self, echr, eloc, udist, ddist):
        try:
            int(echr)
            chr_int_search = True
        except ValueError:
            chr_int_search = False

        if chr_int_search:
            regexstr = "[a-z_]{0,}" + echr
            self.genes.extend(self.session.query(Gene).filter(Gene.seqname.op('regexp')(regexstr))
                              .filter(Gene.strand == '+')
                              .filter(Gene.start.between(eloc - udist, eloc + ddist) |
                                      Gene.end.between(eloc - udist, eloc + ddist) |
                                      ((Gene.start <= eloc) & (Gene.end >= eloc))))
            self.genes.extend(self.session.query(Gene).filter(Gene.seqname.op('regexp')(regexstr))
                              .filter(Gene.strand == '-')
                              .filter(Gene.start.between(eloc - ddist, eloc + udist) |
                                      Gene.end.between(eloc - ddist, eloc + udist) |
                                      ((Gene.start <= eloc) & (Gene.end >= eloc))))
        else:
            self.genes.extend(self.session.query(Gene).filter(Gene.seqname == echr)
                              .filter(Gene.strand == '+')
                              .filter(Gene.start.between(eloc-udist, eloc+ddist) |
                                      Gene.end.between(eloc-udist, eloc+ddist) |
                                      ((Gene.start <= eloc) & (Gene.end >= eloc))))
            self.genes.extend(self.session.query(Gene).filter(Gene.seqname == echr)
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

    def write_results(self, outpath=None, depth=0):
        if outpath is not None:
            ostream = open(outpath, 'w')
            isfile = True
        else:
            ostream = sys.stdout
            isfile = False

        ostream.write("Chromosome\tLoc_start\tLoc_end\tOrientation\tAGI\tType\tAttributes\n")
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
