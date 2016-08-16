# -*- coding: utf-8 -*-

"""tairdbsuite.tairdb_models: module defining the classes which are mapped to the database by sqlalchemy."""

from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql.sqltypes import UnicodeText

Base = declarative_base()


class TairGene(Base):
    __tablename__ = 'tair10genes'
    id = Column(Integer, primary_key=True)
    db = Column(String)
    chromosome = Column(String)
    type = Column(String)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    orientation = Column(String(1))
    agi = Column(String)
    shortsym = Column(String)
    longsym = Column(String)
    attributes = Column(String)
    children = relationship("TairLevel1", back_populates="parent")


class TairLevel1(Base):
    __tablename__ = 'tair10level1'

    id = Column(Integer, primary_key=True)
    model = Column(Integer)
    type = Column(String)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    attributes = Column(String)
    parent_id = Column(Integer, ForeignKey('tair10genes.id'))
    parent = relationship("TairGene", back_populates="children")
    children = relationship("TairLevel2")
    description = relationship("TairDesc", uselist=False, back_populates="level1")


class TairLevel2(Base):
    __tablename__ = 'tair10level2'

    id = Column(Integer, primary_key=True)
    type = Column(String)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    attributes = Column(String)
    level1_id = Column(Integer, ForeignKey('tair10level1.id'))
    level1 = relationship("TairLevel1", back_populates="children")


class TairDesc(Base):
    __tablename__ = 'tair10desc'

    id = Column(Integer, primary_key=True)
    shortdesc = Column(UnicodeText)
    curatorsummary = Column(UnicodeText)
    longdesc = Column(UnicodeText)
    type = Column(String)
    model = Column(Integer)
    level1_id = Column(Integer, ForeignKey('tair10level1.id'))
    level1 = relationship("TairLevel1", back_populates="description")
