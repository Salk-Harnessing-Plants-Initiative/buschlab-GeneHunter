# -*- coding: utf-8 -*-

"""tairdbsuite.tairdb_models: module defining the classes which are mapped to the database by sqlalchemy."""

from sqlalchemy import Column, Integer, Float, String, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql.sqltypes import UnicodeText

Base = declarative_base()


class Gene(Base):
    __tablename__ = 'genes'

    internal_id = Column(Integer, primary_key=True)

    seqname = Column(String)
    source = Column(String)
    feature = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    score = Column(String)
    strand = Column(String(1))
    frame = Column(String(1))
    attribute = Column(String)
    id = Column(String)

    rna = relationship("RNA", back_populates="parent")


class RNA(Base):
    __tablename__ = 'rna'

    internal_id = Column(Integer, primary_key=True)

    seqname = Column(String)
    source = Column(String)
    feature = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    score = Column(String)
    strand = Column(String(1))
    frame = Column(String(1))
    attribute = Column(String)
    id = Column(String)

    parent_id = Column(Integer, ForeignKey('genes.internal_id'))
    parent = relationship("Gene", back_populates="rna")
    features = relationship("Feature", back_populates="parent")


class Feature(Base):
    __tablename__ = 'features'

    internal_id = Column(Integer, primary_key=True)

    seqname = Column(String)
    source = Column(String)
    feature = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    score = Column(String)
    strand = Column(String(1))
    frame = Column(String(1))
    attribute = Column(String)
    id = Column(String)

    parent_id = Column(Integer, ForeignKey('rna.internal_id'))
    parent = relationship("RNA", back_populates="features")


# class TairDesc(Base):
#     __tablename__ = 'tair10desc'
#
#     id = Column(Integer, primary_key=True)
#     shortdesc = Column(UnicodeText)
#     curatorsummary = Column(UnicodeText)
#     longdesc = Column(UnicodeText)
#     type = Column(String)
#     model = Column(Integer)
#     level1_id = Column(Integer, ForeignKey('tair10level1.id'))
#     level1 = relationship("Tair_Level1", back_populates="description")
