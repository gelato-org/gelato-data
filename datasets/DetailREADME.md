Datasets

Variables and associated values coded for Gelato "population samples" are organized into datasets, according to their source (e.g. Autosomal, mtDNA, Y chromosome).

Each dataset is independent, and organized in separate folders. 

    variables.csv: The list of variables for each population sample, coded in a dataset; must contain columns
        SamplePopID: GELATO-wide unique identifier for the sample
        populationName: commonly used name for the sample
        samplesize: number of individuals per each sample
        geographicRegion: continental region, with a coherent genetic history
        dataSet.of.origin: referring to different publication sources, or standardized genetic panels
        lat: latitude of the sampling location (might differ from the glottocode location)
        lon: latitude of the sampling location (might differ from the glottocode location)
        location: relative to the original genetic publication
      glottocode: match with an existing glottocode. 
      languoidName: associated with the glottocode
        curation_notes: notes about the glottocode match from the dataset curators
        Exclude: notes about samples that might be excluded for a clean linguistic/genetic analysis: ok: good language/gene match. exclude: no language/gene match, not enough information to give a match. ADMIXED: the sample is characterized by strong admixture component and might result in difficult comparisons with the rest of the dataset. 
        
      
        
        
    data.csv: The calculated genetic values - i.e. summary statistics. must contain columns
        ExpectedHeterozygosity: Expected heterozygosity (Hs), equivalent to Nei’s genetic diversity [D]
