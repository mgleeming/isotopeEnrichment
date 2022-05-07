# isotopeEnrichment

The isotopeEnrichment.py script determines the abundances of infividual peptide isotopes in liquid chraomatography-mass spectrometry data. The resulting data can be used to determine the extent of heavy isotope incorporation into a target molecule in isotope tracer experiments.

**Workflow**

![Workflow](images/IsotopeEnrichment.png)

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

# How it works

isotopeEnrichment.py uses peptide assignment data from the msmsScans.txt file produced by MaxQuant to calculate the exact mass of theoretical isotopologues that could be expected if the peptide were to be enriched with an arbitrary number of 'heavy' labelled atoms. For each peptide, theoretical exact mases are then used to create an extracted ion chromatogram that is a summation of the intensities of each of the target isotopologues.

A gaussian curve is then fitted to this chromatogram and the target isotopologue intensities are taken as the average of the observed intensities in a given number of mass spectra either side of the gaussian peak maxima.

If a peptide is fully labelled (i.e. 100% isotope incorporation), the entire isotope distribution will be shifted to higher m/z which will likely impair assignment by MaxQuant. To account for this case, we only conduct maxquant searches on the 'control' samples that have not been treated with an isotopic label. The observed retention time and theoretical exact target mass data from peptides assigned by MaxQuant in the control sample are then used to identify the likely corresponding peptide in the isotope-treated samples. While the Gaussian fitting procedure provides some tolerance for drifts in retention time, chromatographic reproducibility in this analysis is of great importance.

In some cases, Gaussian curve fitting may be complicated by the presence of multiple peaks in the extracted window or by low abundance or absent peaks. To ensure that these do not lead to inaccurate or unexpected results, limits can be set on the tolerated range of peak full-width-half-maxima (FWHM). Peaks with FWHMs outside this range will be excluded.

# Examples

To conduct an analysis on all peptides from all proteins identified by MaxQuant in the MS_data_file.mzML data:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt

To conduct an analysis on only peptides from proteins containing '60S' in their name identified by MaxQuant in the MS_data_file.mzML data:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt --searchTerm 60S

To conduct an analysis on only peptides from proteins containing '60S', '40S' or '30S' in their name identified by MaxQuant in the MS_data_file.mzML data:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt --searchTerm 60S --searchTerm 40S --searchTerm 30S

A series of figures can be produced that show the EIC trace for the monoisotope peptide peak as well as the distribution of peptide isotopes determined. Drawing these figures takes additional time and so is disabled by default. Adding the --plot flag will activate plotting features:

    python isotopeEnrichment.py --mzmlFile MS_data_file.mzML --proteinGroupsFile proteinGroups.txt --msmsScansFile msmsScans.txt --searchTerm 60S --plot



# isotopeEnrichment help
    usage: isotopeEnrichment.py [-h] [--searchTerm SEARCHTERM] --mzmlFile MZMLFILE
                                [--proteinGroupsFile PROTEINGROUPSFILE] --modificationSpecificPeptidesFile
                                MODIFICATIONSPECIFICPEPTIDESFILE [--minClusterWidth MINCLUSTERWIDTH]
                                [--addSpecialResidues] [--specialResidue SPECIALRESIDUE]
                                [--avgNSpectra AVGNSPECTRA] [--eicWidth EICWIDTH] [--eicLength EICLENGTH]
                                [--isotopeWeight ISOTOPEWEIGHT] [--processTopN PROCESSTOPN] [--plot]
                                [--plotTopN PLOTTOPN] [--outDirName OUTDIRNAME] [--profile]
                                [--fwhmLim FWHMLIM]

    Extract peptide isotope abundances from LCMS data

    optional arguments:
      -h, --help            show this help message and exit
      --searchTerm SEARCHTERM
                            Text search terms used to filter proteins that are to be indcluded in the
                            analysis. If omitted, all peptides will be included. To specify multiple
                            search terms, include multiple argument/value pairs. For exmaple --searchTerm
                            60S --searchTerm 40S --searchTerm 30S
      --mzmlFile MZMLFILE   File path of mzML data files. To specify multiple mzML files, include multiple
                            argument/value pairs. For example --mzmlFile sample1.mzML --mzmlFile
                            sample2.mzML --mzmlFile sample3.mzML
      --proteinGroupsFile PROTEINGROUPSFILE
                            File path of proteinGroups.txt file produced by MaxQuant
      --modificationSpecificPeptidesFile MODIFICATIONSPECIFICPEPTIDESFILE
                            File path of modificationSpecificPeptides.txt file produced by MaxQuant
      --minClusterWidth MINCLUSTERWIDTH
                            Number of isotopologues to consider in the analysis
      --addSpecialResidues  For each peptide, add the count of residues specified by the --specialResidue
                            flag to the target cluster width.
      --specialResidue SPECIALRESIDUE
                            Creates a column in the output tables containing the total number of
                            occurrences of the specified residue. To specify multiple residues, include
                            multiple argument/value pairs. For example "--specialResidue S
                            --specialResidue G" will create a column with the number of Gly and Ser
                            residues
      --avgNSpectra AVGNSPECTRA
                            Number of spectra either side of the EIC peak maximum to average when
                            calculating isotope abundances
      --eicWidth EICWIDTH   width (in m/z) used to produce EIC plots
      --eicLength EICLENGTH
                            Time range (in min) surrounding a target to produce EIC plots
      --isotopeWeight ISOTOPEWEIGHT
                            mass increment of isotope of interest
      --processTopN PROCESSTOPN
                            Process only the top N most abundant target peptides (measured by MS1
                            precursor intensity).
      --plot                Draw peptide EIC and isotope intensity graphics. This is very time consuming
                            if the number of target peptides is large.
      --plotTopN PLOTTOPN   Draw peptide EIC and isotope intensity graphics for only the top N most
                            abundant target peptides (measured by MS1 precursor intensity). Only active if
                            the --plot flag is given
      --outDirName OUTDIRNAME
                            Results directory name
      --profile             Calculate peak fitting statistics and produce summary plots
      --fwhmLim FWHMLIM     If specified, peptides with a FWHM greater than this value will be ignored.


# Enrichment Calculations (isotopeEnrichment.R & IsoCorrectoR::IsoCorrection())

**A function to optimize peptide isotopolog intensities and transform them into IsoCorrectoR input files**

The isotopeEnrichment.R function arranges peptides into the required input files for natural isotopic abundance (NIA) correction and subsequent mean enrichment calculation. The calculations are performed using the IsoCorrectoR R package:

*Heinrich, P., Kohler, C., Ellmann, L., Kuerner, P., Spang, R., Oefner, P. J., and Dettmer, K. (2018). Correcting for natural isotope abundance and tracer impurity in MS-, MS/MS- and high-resolution-multiple-tracer-data from stable isotope labeling experiments with IsoCorrectoR. Sci. Rep. 8.*

This function allows users to grab the output "intensities.tsv" table from isotopeEnrichment.py and transform it into the correct input format for IsoCorrectoR while selecting and optimal number of isotopolog peaks per peptide entry. Optimization of the isotopolog number relies on estimating first how many isotopologs can be expected from the atomic composition of each peptide and their natural isotopic abundance and secondly from the labelling percentage in soluble amino acid pools and the number of labelled amino acid residues in each peptide sequence.



    ## preparing results files for correction

    IsEnr <- isotopeEnrichment(PyResultsDir = "intensities.tsv",
                               returnCSV = T,
                               verbose = T,
                               rmIsotopologs = 0,
                               OptimizeIsotopologNr = T,
                               files2correct = list.files(path = ".", pattern = "Norm_Factor"),
                               AA4correction = c("Serine", "Glycine"),
                               AAinterprtFileDir = "AminoAcidNames2SingleLetters.csv",
                               ElementalFileDir = "ElementFile.csv",
                               ProtPTMs = c("OX", "AC"),
                               LabelledSamplesNr = 6)

    ## performing NIA correction and mean enrichment calculation

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

# Statistical Filters (EnrichmentSet.R)

**A function to apply thresholds and test statistically peptide enrichment in order to obtain good quality labelled peptide subsets**

The statistical filters applied in this section are meant to remove falsely interpreted isotopolog abundances derived from "labelled" controls, this phenomenon can result from peptide coelution and contamination if the heavy isotopolog peaks with different peptides. This in turn results in an increased relative isotope abundance that does not come from the labelling experiment. This special scenario exemplifies the utility of having non-labelled controls in your samples. This function allows users to apply thresholds on residual labelling, noise and multiple parameters that leverage on the quality of the information that is delivered. The function returns a list of peptide subsets that fit the selected criteria and may be fed as input to the subsequent function in order to annotate their parent protein identities in the experimental dataset.

The filters are a dependendency of the R package [RandoDiStats](https://github.com/MSeidelFed/RandodiStats_package) and their usage and documentation can be found in the provided link.


    ## data reduction based on class comparison algorithms
    
    ## The enrichment file directory must be updated with the new date when run...

    Subset <- EnrichmentSet(EnrichmentFileDir = "IsoCorrectoR_result_MeanEnrichment.csv",
                                      Treatment = as.factor(c(rep("Control_NL",3),
                                                              rep("Control_L",3),
                                                              rep("Treatment_NL",3),
                                                              rep("Treatment_L",3))),
                                      LabelFactor = as.factor(rep(c(rep("Control", 3),
                                                                    rep("Labelled", 3)),2)),
                                      NLMAX = 0.002,
                                      LEnrLack = 0.02,
                                      SigLab = 0.05,
                                      SigLabControl = 0.05,
                                      SigLabTreatment = 0.05, 
                                      CorrectLab = T)

The function is interactive so you will need to type in the R console while running it. The resulting object contains significantly enriched peptides.

# Protein identification

**A function to identify and compile protein enrichment information from a fed peptide subset**

This function allows user identify the parent proteins from peptide subsets, highlight the coverage of the peptides in the protein sequence, retrieve and visualize the noise-corrected enrichments (non-corrected LPFs) and return an output table that is necessary for the next function in the workflow, which corrects LPFs.

ID the parent proteins from the labelled / important peptides

    AnnotateProt <- AnnotateProteins(PeptideVector = Reduce(f = union, x = Subset),
                                     Treatment = as.factor(c(rep("Control_NL",3),
                                                             rep("Control_L",3),
                                                             rep("Treatment_NL",3),
                                                             rep("Treatment_L",3))),
                                     LabelFactor = as.factor(rep(c(rep("Control", 3),
                                                             rep("Labelled", 3)),2)),
                                     FileName = "NL2_in_14N_controls",
                                     Path2FASTA = "160517_Hv_IBSC_PGSB_r1_proteins_HighConf_REPR_annotation.fasta",
                                     Path2MQev = "evidence.txt",
                                     EnrichmentFileDir = "IsoCorrectoR_result_MeanEnrichment.csv",
                                     cexSeq = 0.85,
                                     ProtPTMs = c("OX", "AC"),
                                     SeqCharLength = 70,
                                     GroupPeptides = F,
                                     verbose = T,
                                     CorrectLab = T)
									 
# Labelled peptide fraction (LPFs) correction

**A function to turn non-corrected into corrected labelled peptide/protein fractions (Corr LPFs)**

This function allows users to input their non-corrected LPF matrix from "AnnotateProteins.R" and correct the fractional enrichment in individual peptides using the enrichment percentages in amino acid soluble pools from paired treatments. The function returns the corrected matrix, which can be directly used to calculate fractional synthesis rates.

	Correction_LPF <- LPFcorrection(InputNonCorrMat = AnnotateProt,
                                    files2correct = list.files(path = ".",
                                                               pattern = "Mean_Norm_Factor",
                                                               all.files = T,
                                                               full.names = T,
                                                               recursive = F, include.dirs = T),
                                    AA4correction = c("Serine", "Glycine"),
                                    AAinterprtFileDir = "AminoAcidNames2SingleLetters.csv",
                                    ProtPTMs = c("OX", "AC"),
                                    CorrectMeans = F, 
                                    EnrBoundary = 50, 
                                    GroupPeptides = F)
								
# Protein fractional synthesis rates (K_s)

The calculation uses corrected labelled peptide fractions (LPFs) and multiplies them by the relative growth rate times 100 in order to turn individual enrichments into individual protein fractional snythesis rates. First the calculation requieres a estimate of the growth rates. In our exemplary test-case growth rates were calculated in dependance to drw weight and protein accumulation dynamics, the RGR table look like this:

 | Treatment | corr_mean_RGR |           
 -----------| ----------|------------   |
  4°C       | Cold      | 0.018695011   |
  20°C      | Control   | 0.002082894   |

Subsequently the calculation can be simply done by:

	Ks <- matrix(NA, nrow = nrow(Correction_LPF), ncol = 0)

	for (i in 1:nrow(corr_mean_RGR)) {
  
		Treatment_col_runner <- as.matrix(Correction_LPF[,grep(pattern = rownames(corr_mean_RGR)[i],
                                                               colnames(Correction_LPF)[1:6])])
  
		Ks <- cbind(Ks, (apply(Treatment_col_runner, 2, as.numeric) * corr_mean_RGR[i,2] * 100))
  
	}

	colnames(Ks) <- paste0("Ks_", c(rep(rownames(corr_mean_RGR)[1], 3), rep(rownames(corr_mean_RGR)[2], 3)), "(% per mg-DW per h)")

	Final_table <- cbind(Ks, Correction_LPF)
	
The synthesis rates are the bound at the beggining of the corrected LPF table and ready to export.

