import json
import decimal
import pathlib
import collections

from cldfbench import Dataset as BaseDataset, CLDFSpec, CLDFWriter
from pycldf.sources import Sources
from csvw import Datatype
from csvw.dsv_dialects import Dialect


class Dataset(BaseDataset):
    dir = pathlib.Path(__file__).parent
    id = "gelato"

    def cldf_specs(self):  # A dataset must declare all CLDF sets it creates.
        return CLDFSpec(
            dir=self.cldf_dir,
            module="StructureDataset",
            data_fnames=dict(
                LanguageTable='populations.csv',
                ParameterTable='variables.csv',
                ValueTable='data.csv',
            )
        )

    def cmd_download(self, args):
        pass

    def schema(self, cldf):
        t = cldf.add_columns(
            'LanguageTable',
            {
                'name': 'geographicRegion',
                'dc:description': "Geographic location of the populations is based on information "
                                  "on the genetic samples, and not on linguistic information."
            },
            'country',
            {
                "name": 'samplesize',
                "datatype": "integer",
            },
            {
                "name": 'Average_SNP_count',
                "datatype": "float",
            },
            {
                'name': 'LanguageFamily_Glottocode',
                'dc:description': "Glottocode of the top-level language grouping associated with "
                                  "the population. Language isolates have their own glottocode "
                                  "in this column as well."
            },
            {
                'name': 'LanguageFamily',
                'dc:description': "Name of the top-level language grouping associated with "
                                  "the population."
            },
            'curation_notes_linguistics',
            'curation_notes_genetics',
            {
                "name": "Source",
                "propertyUrl": "http://cldf.clld.org/v1.0/terms.rdf#source",
                "datatype": {"base": "string"},
                "separator": ";"
            }
        )
        cldf['LanguageTable'].common_props['dc:description'] = \
            "Rows in this table represent genetic populations mapped to a language. These " \
            "populations constitute the primary unit of investigation in GeLaTo."
        cldf['LanguageTable', 'Glottocode'].common_props['dc:description'] = \
            "Glottocode identifier, which corresponds to the main language spoken by the " \
            "population. This information is recovered from the original genetic publication, " \
            "and it is extrapolated either from direct sampling observation, cultural/linguistic " \
            "self-identification, or geographical characterization, with the assistance of " \
            "linguists and anthropologists."
        cldf.add_columns(
            'ParameterTable',
            {
                "name": "datatype",
                "datatype": "json",
                "dc:description":
                    "GeLaTo provides parameters (aka variables) of two types. "
                    "Functions of *one* population, where values will be atomic measurements like "
                    "numbers, and functions of a pair of populations, where values are mappings "
                    "of population IDs to atomic values. The latter is distinguished by a value "
                    "of `json` in this column and values must be read as JSON objects. "
                    "The set of values for a variable of "
                    "the latter type can be interpreted as value matrix, i.e. as values for the "
                    "cartesian product of the set of populations."
            },
        )
        lid = cldf['ValueTable', 'languageReference']
        lid.name = 'Population_ID'
        lid.common_props['dc:description'] = "Links a value to a population."
        for fk in cldf['ValueTable'].tableSchema.foreignKeys:
            if fk.columnReference == ['Language_ID']:
                fk.columnReference = ['Population_ID']
        cldf['ValueTable', 'Value'].common_props['dc:description'] = \
            "Either a value with an atomic datatype (like number or string) or a JSON serialized " \
            "mapping of population ID to an atomic value. In the latter case, the corresponding " \
            'variable is a function f of two populations and a value like {"ID1": 5, "ID2" 7} is to ' \
            'be interpreted as f(row[Population_ID], ID1) = 5 and f(row[Population_ID], ID2) = 7.'

    def cmd_makecldf(self, args):
        glangs = args.glottolog.api.cached_languoids
        self.schema(args.writer.cldf)
        args.writer.cldf.sources = Sources.from_file(self.etc_dir / 'sources.bib')
        popname2id = {}
        types = {}
        vc = 0
        for d in self.dir.joinpath('datasets').iterdir():
            if d.is_dir() and d.stem != 'Pemberton_AutosomalSTR':
                for s in d.read_csv('samples.csv', dicts=True):
                    popname2id[s['PopName']] = s['SamplePopID']
                    glang = glangs[s['glottocodeBase']]
                    family = glang.family or glang
                    args.writer.objects['LanguageTable'].append(dict(
                        ID=s['SamplePopID'],
                        Name=s['PopName'],
                        Glottocode=s['glottocodeBase'],
                        Latitude=None if s['lat'] == 'NA' else decimal.Decimal(s['lat'].replace(',', '.')),
                        Longitude=None if s['lon'] == 'NA' else decimal.Decimal(s['lon'].replace(',', '.')),
                        geographicRegion=s['geographicRegion'],
                        Comment=s['Notes_for_the_users'],
                        Source=s['publication'].split('&'),
                        Average_SNP_count=float(s['Average SNP count']),
                        samplesize=int(s['samplesize']),
                        country=s['country'],
                        LanguageFamily_Glottocode=family.id,  # s['glottolog.node1'],
                        LanguageFamily=family.name,  # s['LanguageFamily'],
                        curation_notes_linguistics=s['curation_notes_linguistics'],
                        curation_notes_genetics=s['curation_notes_genetics'],
                    ))
                for row in d.read_csv('variables.csv', dicts=True):
                    if row['type']:
                        types[row['VarID']] = Datatype.fromvalue(row['type'])
                        args.writer.objects['ParameterTable'].append(dict(
                            ID=row['VarID'].replace(' ', '_'),
                            Name=row['Variable name'],
                            Description=row['Description'],
                            datatype=row['type'],
                        ))
                    else:
                        args.log.warning('Skipping untyped variable {}'.format(row['Variable name']))

                for row in d.read_csv('data.csv', dialect=Dialect(lineTerminators=['\r']), dicts=True):
                    #SamplePopID,PopName,
                    # Average SNP count,
                    # medianFST,
                    # medianFSTregion,
                    # MedianFSTAdjustedNeighbors,
                    # numberofNeighborswithin1000km,
                    # harmonicMean,
                    # harmonicMean_5perc,
                    # harmonicMean_95perc
                    for k in row:
                        if k in types:
                            vc += 1
                            args.writer.objects['ValueTable'].append(dict(
                                ID=str(vc),
                                Population_ID=row['SamplePopID'],
                                Parameter_ID=k.replace(' ', '_'),
                                Value=row[k],
                            ))
                pdata = collections.defaultdict(lambda: collections.defaultdict(dict))
                for row in d.read_csv('data_pairwise.csv', dialect=Dialect(lineTerminators=['\r']), dicts=True):
                    for k, v in row.items():
                        if v != 'NA':
                            if k in types:
                                try:
                                    pdata[popname2id[row['Pop2']]][k][popname2id[row['Pop1']]] = float(v)
                                except:
                                    print(k, v)
                                    print(row)
                                    raise
                    #Pop2,Pop1,
                    # FST,
                    # case,

                    # popslistemp,

                    # GEOdist,
                    # FstLinear,

                    # FAMILY,
                    # REGION,

                    # GeneticSplitTime,
                    # GeneticSplitTime_5,
                    # GeneticSplitTime_95
                for popId, data in sorted(pdata.items(), key=lambda i: i[0]):
                    for vid, v in sorted(data.items(), key=lambda i: i[0]):
                        vc += 1
                        args.writer.objects['ValueTable'].append(dict(
                            ID=str(vc),
                            Population_ID=popId,
                            Parameter_ID=vid,
                            Value=json.dumps(v),
                        ))
