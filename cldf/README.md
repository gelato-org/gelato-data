<a name="ds-structuredatasetmetadatajson"> </a>

# StructureDataset GEnes and LAnguages TOgether

**CLDF Metadata**: [StructureDataset-metadata.json](./StructureDataset-metadata.json)

**Sources**: [sources.bib](./sources.bib)

A resource for multidisciplinary studies on human genetic and linguistic variation

property | value
 --- | ---
[dc:bibliographicCitation](http://purl.org/dc/terms/bibliographicCitation) | Barbieri et al. 2022. A global analysis of matches and mismatches between human genetic and linguistic histories. PNAS. DOI: 10.1073/pnas.2122084119
[dc:conformsTo](http://purl.org/dc/terms/conformsTo) | [CLDF StructureDataset](http://cldf.clld.org/v1.0/terms.rdf#StructureDataset)
[dc:identifier](http://purl.org/dc/terms/identifier) | https://gelato.clld.org
[dc:license](http://purl.org/dc/terms/license) | https://creativecommons.org/licenses/by-nc/4.0/
[dcat:accessURL](http://www.w3.org/ns/dcat#accessURL) | https://github.com/gelato-org/gelato-data
[prov:wasDerivedFrom](http://www.w3.org/ns/prov#wasDerivedFrom) | <ol><li><a href="https://github.com/gelato-org/gelato-data/tree/947f1c7">gelato-org/gelato-data 947f1c7</a></li><li><a href="https://github.com/glottolog/glottolog/tree/v4.6">Glottolog v4.6</a></li></ol>
[prov:wasGeneratedBy](http://www.w3.org/ns/prov#wasGeneratedBy) | <ol><li><strong>python</strong>: 3.8.10</li><li><strong>python-packages</strong>: <a href="./requirements.txt">requirements.txt</a></li></ol>
[rdf:ID](http://www.w3.org/1999/02/22-rdf-syntax-ns#ID) | gelato
[rdf:type](http://www.w3.org/1999/02/22-rdf-syntax-ns#type) | http://www.w3.org/ns/dcat#Distribution


## <a name="table-datacsv"></a>Table [data.csv](./data.csv)

property | value
 --- | ---
[dc:conformsTo](http://purl.org/dc/terms/conformsTo) | [CLDF ValueTable](http://cldf.clld.org/v1.0/terms.rdf#ValueTable)
[dc:extent](http://purl.org/dc/terms/extent) | 4791


### Columns

Name/Property | Datatype | Description
 --- | --- | --- 
[ID](http://cldf.clld.org/v1.0/terms.rdf#id) | `string` | Primary key
[Population_ID](http://cldf.clld.org/v1.0/terms.rdf#languageReference) | `string` | Links a value to a population.<br>References [populations.csv::ID](#table-populationscsv)
[Parameter_ID](http://cldf.clld.org/v1.0/terms.rdf#parameterReference) | `string` | References [variables.csv::ID](#table-variablescsv)
[Value](http://cldf.clld.org/v1.0/terms.rdf#value) | `string` | Either a value with an atomic datatype (like number or string) or a JSON serialized mapping of population ID to an atomic value. In the latter case, the corresponding variable is a function f of two populations and a value like {"ID1": 5, "ID2" 7} is to be interpreted as f(row[Population_ID], ID1) = 5 and f(row[Population_ID], ID2) = 7.
[Code_ID](http://cldf.clld.org/v1.0/terms.rdf#codeReference) | `string` | 
[Comment](http://cldf.clld.org/v1.0/terms.rdf#comment) | `string` | 
[Source](http://cldf.clld.org/v1.0/terms.rdf#source) | list of `string` (separated by `;`) | References [sources.bib::BibTeX-key](./sources.bib)

## <a name="table-populationscsv"></a>Table [populations.csv](./populations.csv)

Rows in this table represent genetic populations mapped to a language. These populations constitute the primary unit of investigation in GeLaTo.

property | value
 --- | ---
[dc:conformsTo](http://purl.org/dc/terms/conformsTo) | [CLDF LanguageTable](http://cldf.clld.org/v1.0/terms.rdf#LanguageTable)
[dc:extent](http://purl.org/dc/terms/extent) | 397


### Columns

Name/Property | Datatype | Description
 --- | --- | --- 
[ID](http://cldf.clld.org/v1.0/terms.rdf#id) | `string` | Primary key
[Name](http://cldf.clld.org/v1.0/terms.rdf#name) | `string` | 
[Macroarea](http://cldf.clld.org/v1.0/terms.rdf#macroarea) | `string` | 
[Latitude](http://cldf.clld.org/v1.0/terms.rdf#latitude) | `decimal` | 
[Longitude](http://cldf.clld.org/v1.0/terms.rdf#longitude) | `decimal` | 
[Glottocode](http://cldf.clld.org/v1.0/terms.rdf#glottocode) | `string` | Glottocode identifier, which corresponds to the main language spoken by the population. This information is recovered from the original genetic publication, and it is extrapolated either from direct sampling observation, cultural/linguistic self-identification, or geographical characterization, with the assistance of linguists and anthropologists.
[ISO639P3code](http://cldf.clld.org/v1.0/terms.rdf#iso639P3code) | `string` | 
`Language_Name` | `string` | 
`geographicRegion` | `string` | Geographic location of the populations is based on information on the genetic samples, and not on linguistic information.
`country` | `string` | 
`samplesize` | `integer` | 
`Average_SNP_count` | `float` | 
`LanguageFamily_Glottocode` | `string` | Glottocode of the top-level language grouping associated with the population. Language isolates have their own glottocode in this column as well.
`LanguageFamily` | `string` | Name of the top-level language grouping associated with the population.
`curation_notes_linguistics` | `string` | 
`curation_notes_genetics` | `string` | 
[Source](http://cldf.clld.org/v1.0/terms.rdf#source) | list of `string` (separated by `;`) | References [sources.bib::BibTeX-key](./sources.bib)
[Panel_ID](http://cldf.clld.org/v1.0/terms.rdf#contributionReference) | `string` | Populations are defined by a set of samples from a panel.<br>References [panels.csv::ID](#table-panelscsv)

## <a name="table-variablescsv"></a>Table [variables.csv](./variables.csv)

property | value
 --- | ---
[dc:conformsTo](http://purl.org/dc/terms/conformsTo) | [CLDF ParameterTable](http://cldf.clld.org/v1.0/terms.rdf#ParameterTable)
[dc:extent](http://purl.org/dc/terms/extent) | 14


### Columns

Name/Property | Datatype | Description
 --- | --- | --- 
[ID](http://cldf.clld.org/v1.0/terms.rdf#id) | `string` | Primary key
[Name](http://cldf.clld.org/v1.0/terms.rdf#name) | `string` | 
[Description](http://cldf.clld.org/v1.0/terms.rdf#description) | `string` | 
`datatype` | `json` | GeLaTo provides parameters (aka variables) of two types. Functions of *one* population, where values will be atomic measurements like numbers, and functions of a pair of populations, where values are mappings of population IDs to atomic values. The latter is distinguished by a value of `json` in this column and values must be read as JSON objects. The set of values for a variable of the latter type can be interpreted as value matrix, i.e. as values for the cartesian product of the set of populations.
[Panel_ID](http://cldf.clld.org/v1.0/terms.rdf#contributionReference) | `string` | Variables are defined for the set of populations from a panel.<br>References [panels.csv::ID](#table-panelscsv)

## <a name="table-panelscsv"></a>Table [panels.csv](./panels.csv)

property | value
 --- | ---
[dc:conformsTo](http://purl.org/dc/terms/conformsTo) | [CLDF ContributionTable](http://cldf.clld.org/v1.0/terms.rdf#ContributionTable)
[dc:extent](http://purl.org/dc/terms/extent) | 1


### Columns

Name/Property | Datatype | Description
 --- | --- | --- 
[ID](http://cldf.clld.org/v1.0/terms.rdf#id) | `string` | Primary key
[Name](http://cldf.clld.org/v1.0/terms.rdf#name) | `string` | 
[Description](http://cldf.clld.org/v1.0/terms.rdf#description) | `string` | 
[Contributor](http://cldf.clld.org/v1.0/terms.rdf#contributor) | `string` | 
[Citation](http://cldf.clld.org/v1.0/terms.rdf#citation) | `string` | 

