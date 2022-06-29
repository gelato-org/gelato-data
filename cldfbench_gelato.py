import json
import decimal
import pathlib
import collections

from cldfbench import Dataset as BaseDataset, CLDFSpec
from pycldf.sources import Sources
from csvw import Datatype
from csvw.dsv_dialects import Dialect


class Dataset(BaseDataset):
    dir = pathlib.Path(__file__).parent
    id = "gelato"

    def cldf_specs(self):  # A dataset must declare all CLDF sets it creates.
        return CLDFSpec(dir=self.cldf_dir, module="StructureDataset")

    def cmd_download(self, args):
        pass

    def schema(self, cldf):
        t = cldf.add_component(
            'LanguageTable',
            'geographicRegion',
            'country',
            {
                "name": 'samplesize',
                "datatype": "integer",
            },
            {
                "name": 'Average_SNP_count',
                "datatype": "float",
            },
            'glottolog.node1',
            'LanguageFamily',
            'curation_notes_linguistics',
            'curation_notes_genetics',
            {
                "name": "Source",
                "propertyUrl": "http://cldf.clld.org/v1.0/terms.rdf#source",
                "datatype": {"base": "string"},
                "separator": ";"
            }
        )
        t.common_props['dc:description'] = \
            "Rows in this table represent genetic population mapped to a language. These " \
            "populations constitute the primary unit of investigation in GeLaTo."
        cldf.add_component(
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
        cldf['ValueTable', 'Language_ID'].common_props['dc:description'] = \
            "Links a value to a population."

    def cmd_makecldf(self, args):
        self.schema(args.writer.cldf)
        args.writer.cldf.sources = Sources.from_file(self.etc_dir / 'sources.bib')
        popname2id = {}
        types = {}
        vc = 0
        for d in self.dir.joinpath('datasets').iterdir():
            if d.is_dir() and d.stem != 'Pemberton_AutosomalSTR':
                for s in d.read_csv('samples.csv', dicts=True):
                    popname2id[s['PopName']] = s['SamplePopID']
                    args.writer.objects['LanguageTable'].append(dict(
                        ID=s['SamplePopID'],
                        Name=s['PopName'],
                        Glottocode=s['glottocodeBase'],
                        Latitude=None if s['lat'] == 'NA' else decimal.Decimal(s['lat'].replace(',', '.')),
                        Longitude=None if s['lon'] == 'NA' else decimal.Decimal(s['lon'].replace(',', '.')),
                        geographicRegion=s['geographicRegion'],
                        Comment=s['Notes_for_the_users'],
                        Source=s['publication'].split('&'),
                        Average_SNP_Count=float(s['Average SNP count']),
                        samplesize=int(s['samplesize']),
                    ))
                for row in d.read_csv('variables.csv', dicts=True):
                    """
                    Average SNP count,Average SNP count,Number of successfully typed SNP positions averaged over the individuals within a population. This parameter is proportional to the genotyping quality.
medianFST,Median FST,Median of the FST values distribution associated to each population against all the other populations of the dataset (global value)
medianFSTregion,Median FST within macro region,Median of the FST values distribution associated to each population against the other populations of the same macro geographic region
MedianFSTAdjustedNeighbors,Median FST within geographic neighbors,Median of the FST values distribution associated to each population against neighboring populations (minimum 2) within 1000 km radius and excluding drifted populations
numberofNeighborswithin1000km,Number of neighbors within 1000km,Number of neighbors within a 1000km radius
harmonicMean,Ne - Effective population size,Effective population size calculated with software IBDNe corresponding to the harmonic mean of the last 50 generations 
harmonicMean_5perc,Ne - Effective population size 0.05 percentile,Lower 5% confidence interval associated to the effective population size calculations
harmonicMean_95perc,Ne - Effective population size 0.95 percentile,Upper 5% confidence interval associated to the effective population size calculations

## Model the following with datatype JSON, values being JSON objects mapping Language_ID to val!
FST,FST,Genetic distance between a pair of populations. Calculated with the Weir and Cockerham formula as implemented in PLINK software. Script available at https://github.com/epifaniarango/Fst_forLargeDatasets
FstLinear,Linearized FST,Linearized FST formula as FST/(1-FST)
GEOdist,Geographic distance,Geographic distance between population pairs in km. Calculated from population coordinates with rdist.earth package fields in R (https://www.rdocumentation.org/packages/fields/versions/13.3/topics/rdist.earth)
case,Double or single pair combination,"Each population pair combination is present twice - as Pop1&Pop2 or as Pop2&Pop1 (like in a full distance matrix). Select ""single"" for using only the lower diagonal of the distance matrix."
GeneticSplitTime,Genetic divergence time,Approximate divergence time calculated from Ne with the formula 2Ne x  linearized FST 
GeneticSplitTime_5,Genetic divergence time 0.05 percentile,Lower 5% confidence interval associated to the genetic divergence time
GeneticSplitTime_95,Genetic divergence time 0.95 percentile,Upper 5% confidence interval associated to the genetic divergence time
                    """
                    if row['type']:
                        types[row['VarID']] = Datatype.fromvalue(row['type'])
                        args.writer.objects['ParameterTable'].append(dict(
                            ID=row['VarID'].replace(' ', '_'),
                            Name=row['Variable name'],
                            Description=row['Description'],
                            datatype=row['type'],
                        ))

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
                                Language_ID=row['SamplePopID'],
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
                            Language_ID=popId,
                            Parameter_ID=vid,
                            Value=json.dumps(v),
                        ))
