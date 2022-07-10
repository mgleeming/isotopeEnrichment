import os, sys, pymzml, pickle, shutil, argparse, itertools
import numpy as np
import pandas as pd
import scipy.signal
import scipy.optimize as opt
import pyteomics.mass as pymass
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


DEFAULT_MIN_CLUSTER_WIDTH = 5
DEFAULT_EXTRACTION_WIDTH = 10
DEFAULT_EXTRACTION_LENGTH = 1
DEFAULT_ISOTOPE_WEIGHT = 1
DEFAULT_SPECTRA_TO_AVERAGE = 2
DEFAULT_OUT_DIR_NAME = 'results'
DEFAULT_OUT_FILE_NAME = 'results.dat'

parser = argparse.ArgumentParser(
    description = 'Extract peptide isotope abundances from LCMS data'
)

parser.add_argument('--searchTerm',
                    required = False,
                    action = 'append',
                    type = str,
                    help = 'Text search terms used to filter proteins that are to be indcluded in the analysis. \
                            If omitted, all peptides will be included. To specify multiple search terms, include \
                            multiple argument/value pairs. For exmaple --searchTerm 60S --searchTerm 40S --searchTerm 30S'
                    )
parser.add_argument('--mzmlFileDir',
                    type = str,
                    help = 'Path to directory containing mzML files'
                    )
parser.add_argument('--proteinGroupsFile',
                    type = str,
                    help = 'File path of proteinGroups.txt file produced by MaxQuant'
                    )
parser.add_argument('--modificationSpecificPeptidesFile',
                    type = str,
                    help = 'File path of modificationSpecificPeptides.txt file produced by MaxQuant'
                    )
parser.add_argument('--minClusterWidth',
                    default = DEFAULT_MIN_CLUSTER_WIDTH,
                    type = int,
                    help = 'Number of isotopologues to consider in the analysis'
                    )
parser.add_argument('--addSpecialResidues',
                    action = 'store_true',
                    help = 'For each peptide, add the count of residues specified by the --specialResidue flag to the target cluster width.'
                    )
parser.add_argument('--specialResidue',
                    action = 'append',
                    type = str,
                    help = 'Creates a column in the output tables containing the total number of occurrences of the specified residue. To specify multiple residues, include multiple argument/value pairs. For example "--specialResidue S --specialResidue G" will create a column with the number of Gly and Ser residues'
                    )
parser.add_argument('--avgNSpectra',
                    default = DEFAULT_SPECTRA_TO_AVERAGE,
                    type = int,
                    help = 'Number of spectra either side of the EIC peak maximum to average when calculating isotope abundances'
                    )
parser.add_argument('--eicWidth',
                    default = DEFAULT_EXTRACTION_WIDTH,
                    type = float,
                    help = 'width (in ppm) used to produce EIC plots'
                    )
parser.add_argument('--eicLength',
                    default = DEFAULT_EXTRACTION_LENGTH,
                    type = float,
                    help = 'Time range (in min) surrounding a target to produce EIC plots'
                    )
parser.add_argument('--isotopeWeight',
                    default = DEFAULT_ISOTOPE_WEIGHT,
                    type = float,
                    help = 'mass increment of isotope of interest'
                    )
parser.add_argument('--processTopN',
                    type = int,
                    help = 'Process only the top N most abundant target peptides (measured by MS1 precursor intensity).'
                    )
parser.add_argument('--plot',
                    action = 'store_true',
                    help = 'Draw peptide EIC and isotope intensity graphics. This is very time consuming if the number of target peptides is large.'
                    )
parser.add_argument('--plotTopN',
                    type = int,
                    help = 'Draw peptide EIC and isotope intensity graphics for only the top N most abundant target peptides (measured by MS1 precursor intensity). Only active if the --plot flag is given'
                    )
parser.add_argument('--outDirName',
                    help = 'Results directory name',
                    default = DEFAULT_OUT_DIR_NAME
                    )
parser.add_argument('--outFileName',
                    help = 'Results file name',
                    default = DEFAULT_OUT_FILE_NAME
                    )
parser.add_argument('--profile',
                    action = 'store_true',
                    help = 'Calculate peak fitting statistics and produce summary plots',
                    )
parser.add_argument('--fwhmLim',
                    type = float,
                    help = 'If specified, peptides with a FWHM greater than this value will be ignored.'
                    )
parser.add_argument('--isSpectronaut',
                    action = 'store_true',
                    help = 'Input file is from spectronaut',
                    )
parser.add_argument('--spectronaut_file',
                    type = str,
                    help = 'File path for spectronaut output.'
                    )
parser.add_argument('--writeSchema',
                    action = 'store_true',
                    help = 'Write schema for input file list. Can be used to define treatment and control sample groupings'
                    )
parser.add_argument('--schema',
                    type = str,
                    help = 'Path to the schema file'
                    )

PROTON = 1.00727647

class Quantities(object):

    def __init__(self, mzmlFile, peptide):
        self.mzmlFile = mzmlFile
        self.rts = []
        self.spectrum_indicies = []
        self.clusterWidth = peptide.clusterWidth
        self.intDict = {c:[] for c in range(self.clusterWidth)}

        self.scans_to_extract = None
        self.local_scans_to_extract = None
        return

    def get_score(self):
        try:
            return self.correlation_score_fraction
        except:
            return 0

    def get_rt_boundaries_for_integrated_peak(self):

        try:
            lbound = self.rts[min(self.local_scans_to_extract)]
            rbound = self.rts[max(self.local_scans_to_extract)]
        except TypeError:
            # occurs when scans to extract is None
            return None, None
        return lbound, rbound

    def get_quantities_to_report(self):
        result = []

        if len(self.scans_to_extract) == 0:
            self.result =  [0 for _ in range(self.clusterWidth)]
            return self.result, 0

        for c in range(self.clusterWidth):
            result.append(
                sum(
                    self.intDict[c][ min(self.local_scans_to_extract) : max(self.local_scans_to_extract) ]
                )
            )

        n_spectra_summed = len(self.local_scans_to_extract)
        self.result = result
        return result, n_spectra_summed

    def append_spectrum_index(self, index):
        self.spectrum_indicies.append(index)

    def append_rt(self, rt):
        self.rts.append(rt)

    def append_intensity(self, isotope, intensity):
        self.intDict[isotope].append(intensity)

    def find_peak(self):

        TO_USE = [0,1,2]
        WINDOW = 3
        THRESHOLD = 0.7

        corr_indicies, corr_rts, corr_corrs = [], [], []
        for i in range(len(self.rts)):

            subset = self.rts[i-WINDOW:i+WINDOW]
            if len(subset) < WINDOW:
                corr_rts.append(self.rts[i])
                corr_corrs.append(0)
                #corr_indicies.append(self.spectrum_indicies[i])
                corr_indicies.append(i)
                continue

            traces = np.array([self.intDict[trace][i-WINDOW:i+WINDOW] for trace in TO_USE])

            # pearsons correlation
            corr = np.corrcoef(traces)

            # diagonals are correlation of x with x --> = 1
            # exclude these
            np.fill_diagonal(corr, 0)

            # get total correlation score
            summed_correlation = np.sum(corr)
            if np.isnan(summed_correlation):
                summed_correlation = 0

            corr_rts.append(self.rts[i])
            corr_corrs.append(summed_correlation)
#            corr_indicies.append(self.spectrum_indicies[i])
            corr_indicies.append(i)

        # threshold for total score
        threshold = max(corr_corrs) * 0.5

        if max(corr_corrs) <= 0:
            self.scans_to_extract = []
            return

        corr_rts = np.asarray(corr_rts)
        corr_corrs = np.asarray(corr_corrs)
        corr_indicies = np.asarray(corr_indicies)

        mask = np.where(corr_corrs > threshold)

        # mzml file spectrum indicies of scans above correlation threshold
        corr_scans_threshold = sorted(list(corr_indicies[mask]))

        # the above may contain outliers
        # when a real peak is found, there should be a run of spectra with score > threshold
        # find largest set of consecutive spectra with score > threshold
        scores = []
        for i in range(len(corr_scans_threshold)):
            i_score = 0
            for j in range(i, len(corr_scans_threshold)):
                if (corr_scans_threshold[j] - corr_scans_threshold[i]) == j-i:
                    i_score += 1
            scores.append([i, i_score])

        scores.sort(key=lambda x: x[1], reverse = True)

        # index == index of first spectrum in set
        # step == offset betwen first and last spectrum in set
        index, step = scores[0]

        # scans to extract are scan indicies in the original mzML file
        self.scans_to_extract = np.asarray(np.asarray(self.spectrum_indicies)[mask][index: index+step])

        # local scans to extract are used to find subset of extracted points used for quant
        self.local_scans_to_extract = np.asarray(corr_indicies[mask][index: index+step])

        # find total score for identified region
        self.correlation_score = np.sum(corr_corrs[mask][index: index+step]) / len(corr_corrs[mask][index: index+step])

        max_score = len(TO_USE) * len(TO_USE) - len(TO_USE)
        self.correlation_score_fraction = self.correlation_score / float(max_score)

        # for plotting
        plot = False
        if plot:
            fig = plt.figure()
            axs = fig.subplots(2,1, gridspec_kw={'height_ratios': [3, 1]}, sharex = True)

            for c in range(self.clusterWidth):
                if c in TO_USE:
                    axs[0].plot(self.rts, self.intDict[c])

            min_subset_index = self.spectrum_indicies.index(min(self.scans_to_extract))
            max_subset_index = self.spectrum_indicies.index(max(self.scans_to_extract))
            axs[1].axvspan(self.rts[min_subset_index], self.rts[max_subset_index], color='blue', alpha=0.2)
            axs[1].plot(corr_rts, corr_corrs)
            plt.show()
        return

class Peptide(object):

    def __init__(self, row, options, all_samples, peptide_id):

        self.peptide_id = peptide_id

        if options.isSpectronaut:
            self.read_spectronaut(row, options)
        else:
            self.read_mq(row)

        self.getFormula()
        self.getSpecialResidueCount(options)
        self.getTargets(options)

        self.all_mzml_files = all_samples

        self.peptideQuant = {}
        for mzmlFile in all_samples:
            self.peptideQuant[mzmlFile] = Quantities(mzmlFile, self)

        return

    def get_quantities_to_report(self, mzml_file):
        return self.peptideQuant[mzml_file].get_quantities_to_report()

    def read_mq(self, row):
        self.rt = float(row['Retention time'])
        self.mz = (float(row['Mass']) + float(row['Charge']) * 1.00727647) / float(row['Charge'])
        self.z = int(row['Charge'])
        self.mods = row['Modifications']
        self.sequence = row['Sequence'].upper()
        self.proteinGroup = row['Proteins']
        self.TIC = int(row['Intensity'])
        return

    def read_spectronaut(self, row, options):
        # RT
        self.rt = float(row['EG.ApexRT'].mean())
        self.rt_obs_ll = float(row['EG.StartRT'].min())
        self.rt_obs_hl = float(row['EG.EndRT'].max())

        self.rt_extraction_ll = self.rt_obs_ll - options.eicLength
        self.rt_extraction_hl = self.rt_obs_hl + options.eicLength

        # Mass and m/z
        self.mass = float(row['FG.PrecMz'].mean()) * float(row['FG.Charge'].mean()) - float(row['FG.Charge'].mean()) * PROTON
        self.mz = float(row['FG.PrecMz'].mean())
        self.z = int(row['FG.Charge'].mean())

        # intensities
        self.TIC = int(row['FG.Quantity'].mean())

        # sequence and mods
        row = row.iloc[0]
        self.mods = row['EG.ModifiedSequence']
        self.sequence = row['EG.StrippedSequence'].upper()
        self.proteinGroup = row['PG.ProteinAccessions']
        return

    def appendQuantity(self, mzmlFile, attribute, quantity):
        self.peptideQuant[mzmlFile][attribute].append(quantity)
        return

    def setQuantity(self, mzmlFile, attribute, quantity):
        self.peptideQuant[mzmlFile][attribute] = quantity
        return

    def getQuantity(self, mzmlFile, attribute):
        return self.peptideQuant[mzmlFile][attribute]

    def getSpecialResidueCount(self, options):

        if options.specialResidue:
            self.specialResidueCount = sum(
                [self.sequence.count(sr.upper()) for sr in options.specialResidue]
            )
        else:
            self.specialResidueCount = 0
        return

    def getTargets(self, options):
        self.targets = []

        self.clusterWidth = options.minClusterWidth

        if options.addSpecialResidues:
            self.clusterWidth += self.specialResidueCount

        for isotope in range(self.clusterWidth):
            isotopeCenter = self.mz + float(isotope) * float(options.isotopeWeight) / float(self.z)
            self.targets.append(
                [
                    isotopeCenter - isotopeCenter / 1000000 * options.eicWidth,
                    isotopeCenter + isotopeCenter / 1000000 * options.eicWidth
                ]
            )
        self.minTargetLL = min([_[0] for _ in self.targets]) - 1
        self.maxTargetHL = max([_[1] for _ in self.targets]) + 1
        return

    def getFormula(self):
        composition = pymass.Composition(self.sequence)
        elements = ['C', 'H', 'N', 'O', 'P', 'S']
        self.formula = ''
        for c in composition:
            self.formula += '%s%s ' %(c, composition[c])
        return

    def find_peaks(self):
        for mzml_fine, quant_object in self.peptideQuant.items():
            quant_object.find_peak()
        return

class InputReader(object):

    def __init__(self, options):

        if not any([options.schema, options.writeSchema]):
            print('No schema options supplied')
            sys.exit()

        if options.writeSchema:
            if options.isSpectronaut:
                self.writeSchemaFromSpectronaut()
            else:
                print('INVALID OPTIONS')
            sys.exit()
        return

    def getPeptides(self, options):

        if options.schema:
            self.readSchemaFile(options)
        else:
            print('No schema file supplied - exiting')
            sys.exit()

        if options.isSpectronaut:
            self.readSpectronaut(options)
        return

    def readSchemaFile(self, options):
        # read control files from schema
        self.schema = pd.read_csv(options.schema)
        return

    def writeSchemaFromSpectronaut(self):

        # read spectronaut dataframe
        df = pd.read_csv(options.spectronaut_file, delimiter='\t')

        files = list(set(df['R.FileName'].tolist()))
        of1 = open('fileSchema.tsv','wt')
        of1.write('File,Group\n'%())
        for f in files:
            of1.write('%s,\n'%f)
        of1.close()
        return

    def readSpectronaut(self, options):

        self.peptides = []

        control_samples = list(set(self.schema[self.schema['Group'] == 'Control']['File'].tolist()))
        all_samples = list(set(self.schema['File'].tolist()))

        # read spectronaut dataframe
        df = pd.read_csv(options.spectronaut_file, delimiter='\t')

        # subset entire df to only files in control samples list
        df = df[df['R.FileName'].isin(control_samples)]

        print('Processing peptide targets')
        counter = 1
        for index, group in df.groupby(['EG.PrecursorId']):

            # quick filter
            if len(group) != len(control_samples): continue
            if len(set(group['R.FileName'].tolist())) != len(control_samples): continue

            p = Peptide( group, options, all_samples, counter)

            self.peptides.append(p)
            counter += 1

        print('Found %s peptides' %len(self.peptides))
        return

    def find_peaks(self):

        self.peptides.sort(key=lambda x: x.TIC, reverse=True)

        for i, p in enumerate(self.peptides):
            if (i % 100) == 0:
                print('Processing peak %s of %s' %(i, len(self.peptides)))
            p.find_peaks()
        return

def extract_isotopologue_EICS(reader, options):

    '''
    New fitting procedure
    ---------------------

    1) extract individual isotope EICs
    2) select a subset of these that will be used in step 3
    3) apply a DIA-like procedure to select RT window
        - the subset of EICs should co-elute - find overlapping regions
    4) select top of peaks from DIA reconstruction

    '''

    '''
        1) Extract individual isotopologue EIC traces
    '''
    files = [_ for _ in os.listdir(options.mzmlFileDir) if '.mzml' in _.lower()]
    files_without_extensions = [_.split('.')[0] for _ in files]

    for mzml_file in reader.peptides[0].all_mzml_files:

        print('Processing mzml file %s' %mzml_file)

        # check file exists
        if mzml_file not in files_without_extensions:
            print('Warning: mzml file %s not found! Skipping...')
            continue

        mzml_file_full_path = os.path.join(options.mzmlFileDir, mzml_file) + '.mzML'
        for spectrumi, s in enumerate(pymzml.run.Reader(mzml_file_full_path)):

            if s.ms_level != 1: continue
            time = s.scan_time_in_minutes()
            mzs = s.mz
            ints = s.i

            for i, p in enumerate(reader.peptides):
                if time < p.rt_extraction_ll: continue
                if time > p.rt_extraction_hl: continue

                p.peptideQuant[mzml_file].append_rt(time)
                p.peptideQuant[mzml_file].append_spectrum_index(spectrumi)

                for isotopei, t in enumerate(p.targets):
                    mzIndices = np.where((mzs > t[0]) & (mzs < t[1]))
                    eic_intensity = np.sum(ints[mzIndices])
                    p.peptideQuant[mzml_file].append_intensity(isotopei, eic_intensity)
    return

def write_report(reader, options):

    output = options.outFileName
    try:
        outPath = os.path.join(os.getcwd(), options.outDirName)
        os.makedirs(outPath)
    except:
        print ('\nOutput directory already exists! Exiting...\n')

    of1 = open(os.path.join(outPath, output),'wt')
    maxClusterWidth = max( [p.clusterWidth for p in reader.peptides] )
    text = '\t'.join(
        ['Peptide ID', 'File', 'Sequence', 'Special Residue Count', 'Formula', 'm/z', 'Retantion Time (min)', 'Charge', 'Modifications', 'Protein Gorup', 'Num Spectra Summed', 'Score'] +
        ['Intensity %s' % str(x) for x in range(maxClusterWidth)]
    )
    of1.write('%s\n'%(text))

    for mzml_file in reader.peptides[0].all_mzml_files:
        for p in reader.peptides:

            quantities, num_spectra_summed = p.get_quantities_to_report(mzml_file)

            score = p.peptideQuant[mzml_file].get_score()

            to_write = [str(x) for x in [p.peptide_id, mzml_file, p.sequence, p.specialResidueCount, p.formula, p.mz, p.rt, p.z, p.mods, p.proteinGroup, num_spectra_summed, score]]
            to_write += [str('%.2f' %_) for _ in quantities]

            text = '\t'.join(to_write)
            of1.write('%s\n'%(text))

    of1.close()
    return

def drawPlots(reader, options):

    figDir = 'figs'

    outPath = os.path.join(os.getcwd(), options.outDirName)
    path = os.path.join(outPath, figDir)
    try:
        os.makedirs(path)
    except:
        pass

    mzmlFiles = reader.peptides[0].all_mzml_files

    for ip, p in enumerate(reader.peptides):

        if options.plotTopN and ip > options.plotTopN:
            break

        fig = Figure(figsize=(8, 12), dpi=200)
        axs = fig.subplots(len(mzmlFiles), 2)

        for axi, mzmlFile in enumerate(sorted(mzmlFiles)):

            if len(mzmlFiles) == 1:
                ax = axs
            else:
                ax = axs[axi]

            name = os.path.basename(mzmlFile)

            # chromatograms
            ax[0].set_title( name )
            ax[0].set_ylabel('Abundance')
            ax[0].set_xlim((
                min(p.peptideQuant[mzmlFile].rts),
                max(p.peptideQuant[mzmlFile].rts),
            ))

            lbound, rbound = p.peptideQuant[mzmlFile].get_rt_boundaries_for_integrated_peak()
            for c in range(p.clusterWidth):
                ax[0].plot(
                    p.peptideQuant[mzmlFile].rts,
                    p.peptideQuant[mzmlFile].intDict[c],
                )
                if all([lbound, rbound]):
                    ax[0].axvline(x = lbound, color = 'red', ls = 'dashed', lw = 0.73)
                    ax[0].axvline(x = rbound, color = 'red', ls = 'dashed', lw = 0.75)

            ax[0].set_yticklabels([])

            # isotope abundances
            quantities = p.peptideQuant[mzmlFile].result
            x = list(range(len(quantities)))

            ax[1].set_title( name )
            ax[1].bar( x, quantities , align = 'center' , alpha = 0.5)

            ax[1].set_ylabel('Abundance')
            ax[1].set_yticklabels([])

            if axi == len(mzmlFiles) - 1:
                ax[0].set_xlabel('Retention Time (min)')
                ax[1].set_xlabel('Isotope')
            else:
                ax[0].set_xticklabels([])
                ax[1].set_xticklabels([])

        fig.tight_layout()
        fig.savefig( os.path.join(path, '%s_%s_rt_%s_mz_%s.png' %(p.peptide_id, p.sequence, p.rt, p.mz)))
    return

def main(options):

    # make output directory
    try:
        outPath = os.path.join(os.getcwd(), options.outDirName)
        os.makedirs(outPath)
    except:
        print ('\nOutput directory already exists! Exiting...\n')
#        sys.exit()

    peptides = InputReader(options)
    peptides.getPeptides(options)
    extract_isotopologue_EICS(peptides, options)

    peptides.find_peaks()
    pickle.dump(peptides, open(os.path.join(outPath, options.outFileName + '.pickle','wb')))

    write_report(peptides, options)

    if options.plot:
        drawPlots(peptides, options)

    return

if __name__ == '__main__':
    options =  parser.parse_args()
    main(options)
