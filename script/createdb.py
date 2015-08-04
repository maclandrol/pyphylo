#!/usr/bin/env python
import argparse
import sqlite3
from Bio import SeqIO
import gzip, os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="createdb - Create a light database from mitochondrial genbank file.")
    parser.add_argument("--genbank", "-G", required=True, dest='genbank', help="Your genbank file")
    parser.add_argument("--sqlitedb", '--database', '-db', dest='database', default="genbank.gb", help="Database name. Only sqlite database supported")
    args = parser.parse_args()
    db_exist = os.path.exists(args.database)

    # support only unzip file right now
    records = SeqIO.index(args.genbank, format="genbank")

    schema = '''
        create table IF NOT EXISTS genbank (
            id integer primary key autoincrement not null,
            accession text not null unique,
            seqid text not null unique,
            gi text,
            complete_genome bit default 0,
            source text,
            organism text,
            lineage text
        );
        '''

    # this will be re-designed if needed
    with sqlite3.connect(args.database) as conn:
        conn.executescript(schema)
        for k, rec in records.items():
            annot = rec.annotations
            desc = rec.description
            conn.execute("insert into genbank (accession, seqid, gi, complete_genome, \
                source, organism, lineage) values (?, ?, ?, ?, ?, ?, ?)", [rec.name, rec.id, annot['gi'], \
                ('complete genome' in desc),annot['source'], annot['organism'], ">".join(annot['taxonomy'])])
        conn.commit()
        print("%s elements inserted in %s"%(conn.total_changes, args.database))
