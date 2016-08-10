# -*- coding: utf-8 -*-

"""tairdbsuite.tairdb_creator: module for database creation."""

import os
import sys
import re
import time
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound

from .TairDBModels import Base, TairGene, TairLevel1, TairLevel2, TairDesc


class TairDBCreator(object):
    def __init__(self, dbpath):
        self.engine = create_engine('sqlite:///' + dbpath, echo=False)
        if not os.path.isfile(dbpath):
            Base.metadata.create_all(self.engine)

        session = sessionmaker()
        session.configure(bind=self.engine)
        self.session = session()

    def create_gff_entries(self, gff_filepath):
        file_lines = sum(1 for _ in open(gff_filepath))
        with open(gff_filepath, 'r') as gff_file:
            linenr = 0
            t0 = time.clock()
            tprev = t0
            sys.stdout.write("Reading gff-file '%s'\n" % gff_filepath)
            for line in gff_file:
                line = line.strip()

                linenr += 1
                tcur = time.clock()
                if tcur - tprev >= 1.0:
                    teta = (tcur - t0) / linenr * (file_lines - linenr)
                    percent_cur = float(linenr) / file_lines * 100.0
                    sys.stdout.write(" %3.2f%% (ETA: %.3fs)\r" % (percent_cur, teta))
                    sys.stdout.flush()
                    tprev = tcur

                cols = line.split('\t')
                if cols[2] == 'chromosome':
                    continue

                (eagi, lvl) = self._extract_agi(cols[8])

                if lvl == 0:
                    gene = TairGene(chromosome=cols[0], db=cols[1], type=cols[2], loc_start=int(cols[3]),
                                    loc_end=int(cols[4]), orientation=cols[6], agi=eagi[0], attributes=cols[8])
                    self.session.add(gene)
                    self.session.commit()

                elif lvl == 1:
                    try:
                        gene = self.session.query(TairGene).filter(TairGene.agi == eagi[0]).one()
                        child = TairLevel1(model=int(eagi[1]), type=cols[2],
                                           loc_start=int(cols[3]), loc_end=int(cols[4]), attributes=cols[8])
                        gene.children.append(child)
                        self.session.commit()
                    # except AttributeError:
                    #                         gene.children=[child]
                    #                         self.session.commit()
                    except MultipleResultsFound:
                        sys.stderr.write("multiple results for %s\n" % eagi[0])
                        continue
                    except NoResultFound:
                        sys.stderr.write("no element found for %s\n" % eagi[0])
                        continue

                elif lvl == 2:
                    try:
                        level1 = self.session.query(TairLevel1).join(TairLevel1.parent) \
                            .filter(TairGene.agi == eagi[0]) \
                            .filter(TairLevel1.model == eagi[1]).one()

                        child = TairLevel2(type=cols[2], loc_start=int(cols[3]), loc_end=int(cols[4]),
                                           attributes=cols[8])

                        level1.children.append(child)

                    except AttributeError:
                        level1.children = [child]
                        self.session.commit()
                    except MultipleResultsFound:
                        sys.stderr.write("multiple results for %s.%s\n" % (eagi[0], eagi[1]))
                        continue
                    except NoResultFound:
                        sys.stderr.write("no element found for %s.%s\n" % (eagi[0].eagi[1]))
                        continue

        sys.stdout.write("\n  ok\n")

    def create_desc_entries(self, desc_filepath):
        file_lines = sum(1 for _ in open(desc_filepath))
        with open(desc_filepath, 'r') as desc_file:
            linenr = 0
            t0 = time.clock()
            tprev = t0
            sys.stdout.write("Reading descriptions-file '%s'\n" % desc_filepath)
            for line in desc_file:
                line = line.strip()

                linenr += 1
                tcur = time.clock()
                if tcur - tprev >= 1.0:
                    teta = (tcur - t0) / linenr * (file_lines - linenr)
                    percent_cur = float(linenr) / file_lines * 100.0
                    sys.stdout.write(" %3.2f%% (ETA: %.3fs)\r" % (percent_cur, teta))
                    sys.stdout.flush()
                    tprev = tcur

                cols = line.split('\t')
                eagi = cols[0].split('.')
                emodel = None
                etype = None
                eshortdesc = None
                ecuratorsum = None
                elongdesc = None
                if len(eagi) == 2:
                    try:
                        emodel = int(eagi[1])
                    except ValueError:
                        emodel = None
                if len(cols) > 1:
                    etype = cols[1]
                if len(cols) > 2:
                    eshortdesc = cols[2].decode('utf-8')
                if len(cols) > 3:
                    ecuratorsum = cols[3].decode('utf-8')
                if len(cols) > 4:
                    elongdesc = cols[4].decode('utf-8')

                try:
                    if emodel is not None:
                        level1 = self.session.query(TairLevel1) \
                            .filter(TairLevel1.parent.has(TairGene.agi == eagi[0])) \
                            .filter(TairLevel1.model == emodel).one()

                        level1.description = TairDesc(model=emodel, type=etype,
                                                      shortdesc=eshortdesc,
                                                      curatorsummary=ecuratorsum,
                                                      longdesc=elongdesc)
                except MultipleResultsFound:
                    sys.stderr.write("multiple results for %s.%s level1\n" % (eagi[0], eagi[1]))
                    continue
                except NoResultFound:
                    sys.stderr.write("no element found for %s.%s level1\n" % (eagi[0], eagi[1]))
                    #                     genes=self.session.query(Tair_Gene).filter(Tair_Gene.agi==eagi[0]).all()
                    #                     lvls=self.session.query(Tair_Level1)\
                    #                     .filter(Tair_Level1.parent.has(Tair_Gene.agi==eagi[0])).all()
                    continue
        self.session.commit()
        sys.stdout.write("\n  ok\n")

    def create_name_entries(self, name_filepath):
        file_lines = sum(1 for _ in open(name_filepath))
        with open(name_filepath, 'r') as name_file:
            linenr = 0
            t0 = time.clock()
            tprev = t0
            sys.stdout.write("Reading gene names '%s'\n" % name_filepath)
            for line in name_file:
                line = line.strip()

                linenr += 1
                tcur = time.clock()
                if tcur - tprev >= 1.0:
                    teta = (tcur - t0) / linenr * (file_lines - linenr)
                    percent_cur = float(linenr) / file_lines * 100.0
                    sys.stdout.write(" %3.2f%% (ETA: %.3fs)\r" % (percent_cur, teta))
                    sys.stdout.flush()
                    tprev = tcur

                cols = line.split('\t')

                eagi = cols[0]
                ssym = None
                lsym = None
                if len(cols) > 1:
                    ssym = cols[1]
                if len(cols) > 2:
                    lsym = cols[2]

                try:
                    gene = self.session.query(TairGene).filter(TairGene.agi == eagi).one()
                    gene.shortsym = ssym
                    gene.longsym = lsym
                    self.session.commit()
                except MultipleResultsFound:
                    sys.stderr.write("multiple results for %s level1\n" % eagi)
                    continue
                except NoResultFound:
                    sys.stderr.write("no element found for %s level1\n" % eagi)
                    continue
        sys.stdout.write("\n ok\n")

    def create_sorf_entries(self, sorf_filepath):
        file_lines = sum(1 for _ in open(sorf_filepath))
        with open(sorf_filepath, 'r') as sorf_file:
            linenr = 0
            t0 = time.clock()
            tprev = t0
            sys.stdout.write("Reading SORF genes '%s'\n" % sorf_filepath)
            for line in sorf_file:
                line = line.strip()

                linenr += 1
                tcur = time.clock()
                if tcur - tprev >= 1.0:
                    teta = (tcur - t0) / linenr * (file_lines - linenr)
                    percent_cur = float(linenr) / file_lines * 100.0
                    sys.stdout.write(" %3.2f%% (ETA: %.3fs)\r" % (percent_cur, teta))
                    sys.stdout.flush()
                    tprev = tcur

                if "sORF" not in line:
                    continue

                cols = line.split('\t')

                echr = cols[1].strip('Ath_chr')

                sorf = TairGene(db='sORF', agi=cols[0], chromosome=echr, type='sorf',
                                loc_start=int(cols[2]), loc_end=int(cols[3]),
                                orientation=cols[4], attributes=cols[5] + ';' + cols[6])
                self.session.add(sorf)
            self.session.commit()
        sys.stdout.write("\n ok\n")

    @staticmethod
    def _extract_agi(istr):
        cols = re.split(';|,', istr)  # istr.split(";")

        lvl = None
        agi = None
        for col in cols:
            if 'Parent=' in col:
                agi = col.strip('Parent=').split('.')
                if len(agi) > 1:
                    lvl = 2
                else:
                    lvl = 1
                break
            elif 'ID=' in col:
                idstr = col.strip('ID=').split('-')
                if len(idstr) > 1:
                    agi = idstr[0].split('.')
                    lvl = 2
                else:
                    agi = idstr[0].split('.')
                    if len(agi) > 1:
                        lvl = 1
                    else:
                        lvl = 0
                break
        return agi, lvl
