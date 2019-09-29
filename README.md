# isotopeEnrichment

The isotopeEnrichment.py script determines the abundances of infividual peptide isotopes in liquid chraomatography-mass spectrometry data. The resulting data can be used to determine the extent of heavy isotope incorporation into a target molecule in isotope tracer experiments.

# Usage

isotopeEnrichment.py requires three different input files to run:

- MaxQuant 'proteinGroups.txt' file
- MaxQuant 'msmsScans.txt' file
- 'mzML' file produced from the raw instrument data file.

The isotopeEnrichment.py script uses the results files produced by MaxQuant to provide peptide and protein targeting information. As such, a MaxQuant search must first be conducted against an appropriate protein database file.

The isotopeEnrichment.py script depends on the following python libraries:

- pandas
- matplotlib
- scipy
- pymzml
- pyteomics


# Examples

To conduct an analysis on all peptides from all proteins identified by MaxQuant in the MS_data_file.mzML data:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt

To conduct an analysis on only peptides from proteins containing '60S' in their name identified by MaxQuant in the MS_data_file.mzML data:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt --searchTerm 60S

To conduct an analysis on only peptides from proteins containing '60S', '40S' or '30S' in their name identified by MaxQuant in the MS_data_file.mzML data:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt --searchTerm 60S --searchTerm 40S --searchTerm 30S

By default, a series of figures are produced that show the EIC trace for the monoisotope peptide peak as well as the distribution of peptide isotopes determined. Drawing these figures takes additional time and can be disabled by adding the --noPlot flag:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt --searchTerm 60S --noPlot



# IsotopeEnrichment help

    usage: isotopeEnrichment.py [-h] [--searchTerm SEARCHTERM] --mzmlFile MZMLFILE
                              --proteinGroupsFile PROTEINGROUPSFILE
                              --msmsScansFile MSMSSCANSFILE
                              [--clusterWidth CLUSTERWIDTH]
                              [--avgNSpectra AVGNSPECTRA] [--eicWidth EICWIDTH]
                              [--isotopeWeight ISOTOPEWEIGHT] [--noPlot]
                              [--outDirName OUTDIRNAME]

    Extract peptide isotope abundances from LCMS data

    optional arguments:
      -h, --help            show this help message and exit
      --searchTerm SEARCHTERM
                            Text search terms used to filter proteins that are to
                            be indcluded in the analysis. If omitted, all peptides
                            will be included. To specify multiple search terms,
                            include multiple argument/value pairs. For exmaple
                            --searchTerm 60S --searchTerm 40S --searchTerm 30S
      --mzmlFile MZMLFILE   File path of mzML data file
      --proteinGroupsFile PROTEINGROUPSFILE
                            File path of proteinGroups.txt file produced by
                            MaxQuant
      --msmsScansFile MSMSSCANSFILE
                            File path of msmsScans.txt file produced by MaxQuant
      --clusterWidth CLUSTERWIDTH
                            Number of isotopologues to consider in the analysis
      --avgNSpectra AVGNSPECTRA
                            Number of spectra either side of the EIC peak maximum
                            to average when calculating isotope abundances
      --eicWidth EICWIDTH   width (in m/z) used to produce EIC plots
      --isotopeWeight ISOTOPEWEIGHT
                            mass increment of isotope of interest
      --noPlot              Don't draw peptide EIC and isotope intensity graphics
      --outDirName OUTDIRNAME
                            Results directory name



