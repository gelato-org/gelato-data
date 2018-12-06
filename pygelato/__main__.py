"""
Main command line interface of the pygelato package.

Like programs such as git, this cli splits its functionality into sub-commands
(see e.g. https://docs.python.org/2/library/argparse.html#sub-commands).
The rationale behind this is that while a lot of different tasks may be
triggered using this cli, most of them require common configuration.

The basic invocation looks like

    gelato [OPTIONS] <command> [args]

"""
from __future__ import unicode_literals, division, print_function
import os
import sys

import openpyxl
from clldutils.clilib import ArgumentParser, ParserError, command
from clldutils.path import Path
from clldutils import licenses
from clldutils.markup import Table
from csvw.dsv import UnicodeWriter

import pygelato
from pygelato.util import data_path
from pygelato.dataset import Dataset


HOME = Path(os.path.expanduser('~'))


class ValidationError(ValueError):
    def __init__(self, msg):
        self.msg = msg
        ValueError.__init__(self, msg)


def is_dataset_dir(d):
    return d.exists() and d.is_dir() \
        and d.name != '_template' and d.joinpath('metadata.json').exists()


def get_dataset(args, name=None):
    name = name or args.args[0]
    dir_ = Path(name)
    if not is_dataset_dir(dir_):
        dir_ = data_path(name, repos=args.gelato_repos)
        if not is_dataset_dir(dir_):
            raise ParserError('invalid dataset spec')
    return Dataset(dir_, glottolog_repos=args.glottolog_repos)


@command()
def download(args):
    """
    Download the raw data for a dataset.

    gelato download DATASET_ID
    """
    get_dataset(args).download()


@command()
def ls(args):
    """
    gelato ls [COLS]+

    column specification:
    - license
    - macroareas
    """
    table = Table('ID', 'Title')
    cols = [col for col in args.args if col in ['license', 'macroareas']]
    tl = 40
    if args.args:
        tl = 25
        table.columns.extend(col.capitalize() for col in cols)
    for d in data_path(repos=Path(args.gelato_repos)).iterdir():
        if is_dataset_dir(d):
            ds = Dataset(d)
            row = [d.name, ds.md['dc:title']]
            for col in cols:
                if col == 'license':
                    lic = licenses.find(ds.md.get('dc:license') or '')
                    row.append(lic.id if lic else ds.md.get('dc:license'))

            table.append(row)
    print(table.render(tablefmt='simple', sortkey=lambda r: r[0], condensed=False))


@command()
def ho(args):
    """
    SamplePopID,populationName,samplesize,geographicRegion,dataSet.of.origin,lat,lon,location,languoidName,glottocode,curation_notes,Exclude

    population	samplesize	geographicRegion	lat	lon	Country	FreqSharingTot	FreqsharingWithin
    FreqSharingOut	Mean_lengthsharingWithin	totallenghtsharing	howmanypops	homozyg	variancehomozyg
    averageFST	averageFSTregion	PANEL	glottocodeBase	glottocodeRecent	glottolog.node1	notes	curationNotes

    """
    m = {
        'SamplePopID': 'SamplePopID',
        'population': 'populationName',
        'samplesize': 'samplesize',
        'geographicRegion': 'geographicRegion',
        #'': 'dataSet.of.origin',
        'lat': 'lat',
        'lon': 'lon',
        'Country': 'location',
        #'': 'languoidName',
        'glottocodeBase': 'glottocode',
        'curationNotes': 'curation_notes',
        #'': 'Exclude',
    }
    m = {v: k for k, v in m.items()}
    sheader = 'SamplePopID populationName samplesize geographicRegion lat lon location glottocode curation_notes'.split()
    sheet = openpyxl.load_workbook(
        str(Path(__file__).parent.parent.parent / 'popmore5individuals_geninfo.xlsx'),
        data_only=True).active
    header = []

    dsdir = Path(__file__).parent.parent / 'datasets' / 'HumanOrigins_AutosomalSNP'
    with UnicodeWriter(dsdir / 'samples.csv') as s, UnicodeWriter(dsdir / 'data.csv', delimiter=';') as d:
        s.writerow(sheader)
        d.writerow('SamplePopID;populationName.x;homozyg;averageFST;averageFSTregion'.split(';'))
        for i, row in enumerate(sheet.rows):
            row = [c.value for c in row]
            if i == 0:
                header = row
                continue
            row = dict(zip(header, row))
            row['SamplePopID'] = str(i)
            s.writerow([row[m[h]] for h in sheader])
            d.writerow([row['SamplePopID'], row['population'], row['homozyg'], row['averageFST'], row['averageFSTregion']])


def main():
    parser = ArgumentParser('pygelato')
    parser.add_argument(
        '--gelato-repos',
        help="path to gelato data repository",
        default=Path(pygelato.__file__).parent.parent)
    parser.add_argument(
        '--glottolog-repos',
        help="path to glottolog data repository",
        default=None)
    sys.exit(parser.main())
