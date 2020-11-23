# isotopeEnrichment

The isotopeEnrichment.py script determines the abundances of infividual peptide isotopes in liquid chraomatography-mass spectrometry data. The resulting data can be used to determine the extent of heavy isotope incorporation into a target molecule in isotope tracer experiments.

# Usage

isotopeEnrichment.py requires three different input files to run:

- MaxQuant 'evidence.txt' file
- 'mzML' file produced from the raw instrument data file.

The isotopeEnrichment.py script uses the results files produced by MaxQuant to provide peptide and protein targeting information. As such, a MaxQuant search must first be conducted against an appropriate protein database file.

The isotopeEnrichment.py script depends on the following python libraries:

- pandas
- matplotlib
- scipy
- pymzml
- pyteomics

# How it works

isotopeEnrichment.py uses peptide assignment data from the evidence.txt file produced by MaxQuant to calculate the exact mass of theoretical isotopologues that could be expected if the peptide were to be enriched with an arbitrary number of 'heavy' labelled atoms. For each peptide, theoretical exact mases are then used to create an extracted ion chromatogram that is a summation of the intensities of each of the target isotopologues.

# Examples

To conduct an analysis on all peptides from all proteins identified by MaxQuant in the MS_data_file.mzML data:

    python isotopeEnrichmentv2.py --mzmlFileDir /path/to/mzml/dir --evidenceFile /path/to/evidence.txt

By default, isotopeEnrichment considers 5 isotopes for any given cluster. If you've used labelled amino acids, and want to increase this limit based on the presece of these:

    python isotopeEnrichmentv2.py --mzmlFileDir /path/to/mzml/dir --evidenceFile /path/to/evidence.txt --specialResidue S --specialResidue G --addSpecialResidues

Sometimes, not all peptides are detected in all samples. To exclude peptides from the analysis that were detected in too few samples, add the --minSamples argument followed by the desired threshold. Setting this parameter equal to the number of mzml files will require every peptide to be present in every sample.

    python isotopeEnrichmentv2.py --mzmlFileDir /path/to/mzml/dir --evidenceFile /path/to/evidence.txt --minSamples 10

A series of figures can be produced that show the EIC trace for the monoisotope peptide peak as well as the distribution of peptide isotopes determined and a representative mass spectrum of the isotope distribution. Drawing these figures takes additional time and so is disabled by default. Adding the --plot flag will activate plotting features. Further adding the --plotTopN flag will plot figures for the top N peptides by total intensity only.

    python isotopeEnrichmentv2.py --mzmlFileDir /path/to/mzml/dir --evidenceFile /path/to/evidence.txt --plot --plotTopN 20




# IsotopeEnrichment help
    usage: isotopeEnrichmentv2.py [-h] --evidenceFile EVIDENCEFILE --mzmlFileDir
                                  MZMLFILEDIR [--minClusterWidth MINCLUSTERWIDTH]
                                  [--addSpecialResidues]
                                  [--specialResidue SPECIALRESIDUE]
                                  [--avgNSpectra AVGNSPECTRA]
                                  [--intWidth INTWIDTH] [--eicLength EICLENGTH]
                                  [--isotopeWeight ISOTOPEWEIGHT] [--plot]
                                  [--plotTopN PLOTTOPN] [--minSamples MINSAMPLES]
                                  [--outDir OUTDIR]

    Extract peptide isotope abundances from LCMS data

    optional arguments:
      -h, --help            show this help message and exit
      --evidenceFile EVIDENCEFILE
                            MaxQuant evidence.txt file
      --mzmlFileDir MZMLFILEDIR
                            Directory containing mzML files. Files in this
                            directory with an .mzML extention (case insensitive)
                            will be matched.
      --minClusterWidth MINCLUSTERWIDTH
                            Number of isotopologues to consider in the analysis
      --addSpecialResidues  For each peptide, add the count of residues specified
                            by the --specialResidue flag to the target cluster
                            width.
      --specialResidue SPECIALRESIDUE
                            Creates a column in the output tables containing the
                            total number of occurrences of the specified residue.
                            To specify multiple residues, include multiple
                            argument/value pairs. For example "--specialResidue S
                            --specialResidue G" will create a column with the
                            number of Gly and Ser residues
      --avgNSpectra AVGNSPECTRA
                            Number of spectra either side of the EIC peak maximum
                            to average when calculating isotope abundances
      --intWidth INTWIDTH   width (in ppm) used to integrate peaks
      --eicLength EICLENGTH
                            Time range (in min) surrounding a target to produce
                            EIC plots
      --isotopeWeight ISOTOPEWEIGHT
                            mass increment of isotope of interest
      --plot                Draw peptide EIC and isotope intensity graphics. This
                            is very time consuming if the number of target
                            peptides is large.
      --plotTopN PLOTTOPN   Draw peptide EIC and isotope intensity graphics for
                            only the top N most abundant target peptides (measured
                            by MS1 precursor intensity). Only active if the --plot
                            flag is given
      --minSamples MINSAMPLES
                            Discard peptides observed in fewer than this number of
                            samples
      --outDir OUTDIR       Results directory name


# enrichment calculations

The isotopeEnrichment.R function arranges peptides into the required input files for natural isotopic abundance (NIA) correction and subsequent mean enrichment calculation. The caluclations are performed using the IsoCorrectoR R package:

*Heinrich, P., Kohler, C., Ellmann, L., Kuerner, P., Spang, R., Oefner, P. J., and Dettmer, K. (2018). Correcting for natural isotope abundance and tracer impurity in MS-, MS/MS- and high-resolution-multiple-tracer-data from stable isotope labeling experiments with IsoCorrectoR. Sci. Rep. 8.*


    ## preparing results files for correction

    IsEnr <- isotopeEnrichment(PyResultsDir = "result.dat", returnCSV = T, verbose = T)

    ## performing NIA correction and mean enrichment calculation

    library(IsoCorrectoR)

    Enrichment <- IsoCorrectoR::IsoCorrection(MeasurementFile = "MeasurementFile.csv",
                                              ElementFile = "ElementFile.csv",
                                              MoleculeFile = "MoleculeFile.csv",
                                              CorrectTracerImpurity = T,
                                              CorrectTracerElementCore = T,
                                              CalculateMeanEnrichment = T,
                                              UltraHighRes = F,
                                              FileOutFormat = "csv",
                                              ReturnResultsObject = T,
                                              CorrectAlsoMonoisotopic = T,
                                              verbose = T)
