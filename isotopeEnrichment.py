import os, sys, pymzml, pickle, shutil, argparse, itertools
import numpy as np
import pandas as pd
import scipy.signal
import scipy.optimize as opt
import pyteomics.mass as pymass
import matplotlib.pyplot as plt

DEFAULT_MIN_CLUSTER_WIDTH = 5
DEFAULT_EXTRACTION_WIDTH = 0.01
DEFAULT_EXTRACTION_LENGTH = 0.5
DEFAULT_ISOTOPE_WEIGHT = 1
DEFAULT_SPECTRA_TO_AVERAGE = 2
DEFAULT_OUT_DIR_NAME = 'results'

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
                    help = 'width (in m/z) used to produce EIC plots'
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
        return

    def append_spectrum_index(self, index):
        self.spectrum_indicies.append(index)

    def append_rt(self, rt):
        self.rts.append(rt)

    def append_intensity(self, isotope, intensity):
        self.intDict[isotope].append(intensity)

    def find_peak(self):
        print("FIND")

        TO_USE = [0,1,2]

        print(len(self.rts))
        print(len(self.spectrum_indicies))
        for c in range(self.clusterWidth):
            print(len(self.intDict[c]))

        for combination in itertools.combinations(TO_USE, 2):
            print(combination)

        # find RMSD for a sliding window

        # set parameter for window length
        # this is the expected peak width in scans

        # find RMSD for window
        # slide window across rt range

        # min RMSD beccomes scan center
        sys.exit()
        return

class Peptide(object):

    def __init__(self, row, options, all_samples ):

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

            # TODO
            # this lists targets as <ISOTOPE> above the ms precursor target
            # -- is this misleading?
            # TODO

            isotopeCenter = self.mz + float(isotope) * float(options.isotopeWeight) / float(self.z)
            self.targets.append(
                [
                    isotopeCenter - options.eicWidth,
                    isotopeCenter + options.eicWidth
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
                self.writeSchemaFromSpectronaut(df)
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

    def writeSchemaFromSpectronaut(self, df):
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
        for index, group in df.groupby(['EG.PrecursorId']):

            # quick filter
            if len(group) != len(control_samples): continue
            if len(set(group['R.FileName'].tolist())) != len(control_samples): continue

            p = Peptide( group, options, all_samples )
            self.peptides.append(p)

        print('Found %s peptides' %len(self.peptides))
        return

    def find_peaks(self):
        for p in self.peptides:
            p.find_peaks()
        return

def gauss(x, amp, cent, wid, scale = 1):
    return(amp/ (np.sqrt(2*np.pi*(wid/scale)**2 )) * np.exp(-(x-(cent/scale))**2 / (2*(wid/scale)**2)))

def getPeptides(target_protein_ids, options):
    peptides = []

    df = pd.read_csv(options.modificationSpecificPeptidesFile, delimiter='\t')

    # filter to only identified spectra
    df = df[df['Reverse'] != '+']
    df = df[df['Potential contaminant'] != '+']

    print('Processing peptide targets')
    for index, row in df.iterrows():

        if target_protein_ids:
            # test if this MSMS spectrum matches a protein in target_protein_ids
            matches = [_ for _ in target_protein_ids if _ in row['Proteins']]
            if len(matches) == 0: continue

        for charge in row['Charges'].split(';'):
            row['Charge'] = charge
            p = Peptide( row, options )
            peptides.append(p)

    peptides.sort(key = lambda x: x.TIC, reverse = True)

    if options.processTopN:
        try:
            return peptides[0:options.processTopN+1]
        except IndexError:
            print('Warning - asked to analyse top %s peptides but only %s were found. Processing all...' %(
                options.processTopN, len(peptides))
            )

    return peptides

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
        break
    return

def getEICData_old(peptides, outPath, options):

    output = 'result.dat'

    print('Found %s peptide targets' %len(peptides))

    of1 = open(os.path.join(outPath, output),'wt')
    maxClusterWidth = max( [p.clusterWidth for p in peptides] )
    text = '\t'.join(
        ['File', 'Sequence', 'Special Residue Count', 'Formula', 'm/z', 'Retantion Time (min)', 'Charge', 'Modifications', 'Protein Gorup'] + ['Intensity %s' % str(x) for x in range(maxClusterWidth)]
    )
    of1.write('%s\n'%(text))

    for mzmlFile in options.mzmlFile:
        print('Processing %s' %mzmlFile)

        '''
        Fitting procedure
        ------------------------------------------------

        1) Create a summed EIC for each target ion
            - ensures that some intensity will exist
              even for completly labeled peptides

        2) Find index of spectrum at EIC peak maximum

        3) For spectra within avgNSpectra of max,
            - average abundance values for each isotope

        '''

        # get summed EIC intensity
        spectra = pymzml.run.Reader(mzmlFile)
        for s in spectra:

            if s.ms_level != 1: continue
            time = s.scan_time_in_minutes()
            mzs = s.mz
            ints = s.i

            for i, p in enumerate(peptides):
                if abs(p.rt - time) < options.eicLength:

                    # used for picking EIC peak max
                    p.appendQuantity(mzmlFile, 'rts', time)

                    summedEIC = 0
                    for t in p.targets:
                        mzIndices = np.where((mzs > t[0]) & (mzs < t[1]))
                        summedEIC += np.sum(ints[mzIndices])
                    p.appendQuantity(mzmlFile, 'ints', summedEIC)

                    # used later for actual isotope abundance calculation

                    # don't need whole spectrum
                    mask = np.where( (mzs > p.minTargetLL) & (mzs < p.maxTargetHL) )
                    mzsubset = mzs[mask]
                    intsubset = ints[mask]

                    p.appendQuantity(mzmlFile, 'mzs', mzsubset)
                    p.appendQuantity(mzmlFile, 'mzInts', intsubset)

        # fit gaussians
        for p in peptides:
            cen = p.rt
            wid = max(p.getQuantity(mzmlFile, 'rts')) - min(p.getQuantity(mzmlFile, 'rts'))
            amp = max(p.getQuantity(mzmlFile, 'ints'))
            p0 = [amp,cen,wid]

            try:
                RTpopt, RTpcov= opt.curve_fit(
                    lambda x, amp, cen, wid: gauss(x, amp, cen, wid, scale = 1),
                    p.getQuantity(mzmlFile, 'rts'), p.getQuantity(mzmlFile, 'ints'), p0=p0, maxfev = 2000
                )
            except RuntimeError:
                # optimal fit parameters not found
                continue

            # calcualte full width at half max
            eicFWHM = 2*np.sqrt(2*np.log(2))*RTpopt[2] / 1

            if options.fwhmLim:
                if eicFWHM > options.fwhmLim:
                    continue

            # save fitting parameters
            p.setQuantity(mzmlFile, 'fits', [RTpopt, RTpcov])
            p.setQuantity(mzmlFile, 'eicFWHM', eicFWHM)

            # use curve to create fitted gaussian intensities at given rts
            fit_ints = gauss(
                p.getQuantity(mzmlFile, 'rts'), *p.getQuantity(mzmlFile, 'fits')[0]
            )

            # pick and save index of maximal intensity
            # -- this is used as the center point for the window of spectra to average
            index = np.argmax(fit_ints)
            p.setQuantity(mzmlFile, 'maxrt', p.getQuantity(mzmlFile, 'rts')[index])
            p.setQuantity(mzmlFile, 'maxrtindex', index)

            abundances = []
            for t in p.targets:
                abundances.append([])
                for index, eicRT in enumerate(p.getQuantity(mzmlFile, 'rts')):
                    if abs(index - p.getQuantity(mzmlFile, 'maxrtindex')) < options.avgNSpectra:
                        specmzs = p.getQuantity(mzmlFile, 'mzs')[index]
                        specints = p.getQuantity(mzmlFile, 'mzInts')[index]
                        mask = np.where(
                            (specmzs > t[0])
                            &
                            (specmzs < t[1])
                        )
                        intsubset = specints[mask]
                        abundances[-1].append(np.sum(intsubset))

            # used only for plotting
            p.setQuantity(mzmlFile, 'abundances', abundances)

            text = '\t'.join(
                [str(x) for x in [mzmlFile, p.sequence, p.specialResidueCount, p.formula, p.mz, p.rt, p.z, p.mods, p.proteinGroup]] + [str(sum(x) / len(x)) for x in p.getQuantity(mzmlFile, 'abundances')]
            )
            of1.write('%s\n'%(text))

    of1.close()
    return

def drawProfilePlots(peptides, outPath, options):

    figDir = 'profile'

    path = os.path.join(outPath, figDir)
    try:
        os.makedirs(path)
    except:
        pass

    plotxlim = 3
    vals = []
    for mzmlFile in options.mzmlFile:
        for p in peptides:
            try:
                v = abs(p.getQuantity(mzmlFile, 'eicFWHM'))
                if v > plotxlim:
                    v = plotxlim
                vals.append(v)
            except:
                pass

    fig, ax = plt.subplots()

    ax.hist(vals, density=False, bins=100, range=(0,plotxlim))
    ax.set_ylabel('Counts')
    ax.set_xlabel('EIC FWHM');

    fig.tight_layout()
    fig.savefig( os.path.join(path, 'FWHM_Histogram.png'))

    return

def drawPlots(peptides, outPath, options):

    figDir = 'figs'

    path = os.path.join(outPath, figDir)
    try:
        os.makedirs(path)
    except:
        pass

    for ip, p in enumerate(peptides):

        if options.plotTopN and ip > options.plotTopN:
            break

        fig, axs = plt.subplots(len(options.mzmlFile), 2)

        for axi, mzmlFile in enumerate(options.mzmlFile):

            if len(options.mzmlFile) == 1:
                ax = axs
            else:
                ax = axs[axi]

            name = os.path.basename(mzmlFile)
            try:
                f = p.getQuantity(mzmlFile, 'fits')[0]
            except KeyError:
                continue

            # chromatograms
            ax[0].set_title( name )
            ax[0].set_ylabel('Abundance')
            ax[0].set_xlim((
                min(p.getQuantity(mzmlFile, 'rts')),
                max(p.getQuantity(mzmlFile, 'rts'))
            ))
            ax[0].plot(
                p.getQuantity(mzmlFile, 'rts'),
                p.getQuantity(mzmlFile, 'ints')
            )

            ax[0].plot(
                p.getQuantity(mzmlFile, 'rts'),
                gauss(p.getQuantity(mzmlFile, 'rts'), *p.getQuantity(mzmlFile, 'fits')[0])
            )
            ax[0].set_yticklabels([])

            # plot chromatogram boundary lines

            HWHM = p.getQuantity(mzmlFile, 'eicFWHM')/2
            lbound = p.getQuantity(mzmlFile, 'maxrt') - HWHM
            rbound = p.getQuantity(mzmlFile, 'maxrt') + HWHM

            ax[0].axvline(x = lbound, color = 'red', ls = 'dashed', lw = 0.73)
            ax[0].axvline(x = rbound, color = 'red', ls = 'dashed', lw = 0.75)
            ax[0].axvline(x = p.rt, color = 'blue', ls = 'dashed', lw = 0.75)

            # isotope abundances
            x = list(range(len(p.getQuantity(mzmlFile, 'abundances'))))
            ax[1].set_title( name )
            ax[1].bar(
                x, [sum(x) for x in p.getQuantity(mzmlFile, 'abundances')] , align = 'center' , alpha = 0.5
            )

            ax[1].set_ylabel('Abundance')
            ax[1].set_yticklabels([])

            if axi == len(options.mzmlFile) - 1:
                ax[0].set_xlabel('Retention Time (min)')
                ax[1].set_xlabel('Isotope')

        fig.tight_layout()
        plt.savefig( os.path.join(path, '%s_%s_rt_%s_mz_%s.png' %(len(p.sequence), p.sequence, p.rt, p.mz)))
        plt.clf()
        plt.cla()
        plt.close()
    return

def main(options):

    # make output directory
    try:
        outPath = os.path.join(os.getcwd(), options.outDirName)
        os.makedirs(outPath)
    except:
        print ('\nOutput directory already exists! Exiting...\n')
#        sys.exit()

#    peptides = getPeptides(options)
    import pickle
    test = True
    if not test:
        peptides = InputReader(options)
        peptides.getPeptides(options)
        extract_isotopologue_EICS(peptides, options)
        pickle.dump(peptides, open('peptides.pickle','wb'))
    else:
        peptides = pickle.load(open('peptides.pickle','rb'))

    peptides.find_peaks()
    sys.exit()

    getEICData(peptides, outPath, options)

    if options.profile:
        drawProfilePlots(peptides, outPath, options)

    if options.plot:
        drawPlots(peptides, outPath, options)

    return

if __name__ == '__main__':
    options =  parser.parse_args()
    main(options)
