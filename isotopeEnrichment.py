import os, sys, pymzml, pickle, shutil, argparse
import numpy as np
import pandas as pd
import scipy.signal
import scipy.optimize as opt
import pyteomics.mass as pymass
import matplotlib.pyplot as plt

DEFAULT_ISOTOPE_COUNTER = 6
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
parser.add_argument('--mzmlFile',
                    required = True,
                    action = 'append',
                    type = str,
                    help = 'File path of mzML data files. To specify multiple mzML files, include multiple \
                            argument/value pairs. For example --mzmlFile sample1.mzML --mzmlFile sample2.mzML \
                            --mzmlFile sample3.mzML'
                    )
parser.add_argument('--proteinGroupsFile',
                    required = True,
                    type = str,
                    help = 'File path of proteinGroups.txt file produced by MaxQuant'
                    )
parser.add_argument('--msmsScansFile',
                    required = True,
                    type = str,
                    help = 'File path of msmsScans.txt file produced by MaxQuant'
                    )
parser.add_argument('--clusterWidth',
                    default = DEFAULT_ISOTOPE_COUNTER,
                    type = int,
                    help = 'Number of isotopologues to consider in the analysis'
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

parser.add_argument('--specialResidue',
                    action = 'append',
                    type = str,
                    help = 'Creates a column in the output tables containing the total number of occurrences of the specified residue. To specify multiple residues, include multiple argument/value pairs. For example "--specialResidue S --specialResidue G" will create a column with the number of Gly and Ser residues'
                    )

class Peptide(object):
    def __init__(self, row, options):

        self.peptideQuant = {}
        for mzmlFile in options.mzmlFile:
            self.peptideQuant[mzmlFile] = {
                'rts' : [], 'ints' : [], 'mzs' : [], 'mzInts' : []
            }

        if pd.isna(row['Precursor apex offset time']):
            offset = 0
        else:
            offset = float(row['Precursor apex offset time'])
        self.rt = float(row['Retention time']) + offset
        self.trigger_rt = float(row['Retention time'])

        self.TIC = float(row['Total ion current'])
        self.mz = float(row['m/z'])
        self.z = int(row['Charge'])
        self.mods = row['Modifications']
        self.mod_seq = row['Modified sequence']
        self.sequence = row['Sequence'].upper()
        self.proteinGroup = row['Proteins']

#        print('')
#        print(self.rt)
#        print(self.trigger_rt)
        self.getTargets(options)
        self.getFormula()
        self.getSpecialResidueCount(options)
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
        for isotope in range(options.clusterWidth):

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

def gauss(x, amp, cent, wid, scale = 1):
    return(amp/ (np.sqrt(2*np.pi*(wid/scale)**2 )) * np.exp(-(x-(cent/scale))**2 / (2*(wid/scale)**2)))

def getPeptides(target_protein_ids, options):
    peptides = []

    df = pd.read_csv(options.msmsScansFile, delimiter='\t')

    # filter to only identified spectra
    df = df[df['Identified'] == '+']

    peptideIndes = {}
    counter = 0

    print('Processing peptide targets')
    for index, row in df.iterrows():

        counter += 1

        # skip unassigned or contaminant peptides
        if len(row['Sequence']) < 2 or len(row['Proteins']) < 2 or 'CON__' in row['Proteins'] or row['Identified'] != '+': continue

        if row['Sequence'] in peptideIndes.keys():
            if abs(peptideIndes[row['Sequence']]['mz'] - float(row['m/z'])) < 0.001:
                if int(peptideIndes[row['Sequence']]['z']) == int(row['Charge']):
                    if abs(peptideIndes[row['Sequence']]['rt'] - float(row['Retention time'])) < 2:
                        continue
        else:
            peptideIndes[row['Sequence']] = {
                'mz': float(row['m/z']),
                'z' : int(row['Charge']),
                'rt' : float(row['Retention time'])
            }

        # test if this MSMS spectrum matches a protein in target_protein_ids
        matches = [_ for _ in target_protein_ids if _ in row['Proteins']]
        if len(matches) == 0: continue

        # skip cases where charge in undetermined
        if int(row['Charge']) == 0: continue

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

def getEICData(peptides, outPath, options):
#    base = os.path.basename(options.mzmlFile)
#    output = 'result_%s.dat' %base.split('.')[0]
    output = 'result.dat'

    print('Found %s peptide targets' %len(peptides))

    of1 = open(os.path.join(outPath, output),'wt')

    text = '\t'.join(
        ['File', 'Sequence', 'Special Residue Count' 'Formula', 'm/z', 'Retantion Time (min)', 'Charge', 'Modifications', 'Protein Gorup'] + ['Intensity %s' % str(x) for x in range(options.clusterWidth)]
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
                [str(x) for x in [mzmlFile, p.specialResidueCount, p.sequence, p.formula, p.mz, p.rt, p.z, p.mods, p.proteinGroup]] + [str(sum(x) / len(x)) for x in p.getQuantity(mzmlFile, 'abundances')]
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
#            try:
#                name = str(p.getQuantity(mzmlFile, 'eicFWHM'))
#            except:
#                name = 'None'
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
            ax[0].axvline(x = p.trigger_rt, color = 'green', ls = 'dashed', lw = 0.75)

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
        sys.exit()

    # get list of proteins containing '60S' string
    df = pd.read_csv(options.proteinGroupsFile, delimiter='\t')

    if options.searchTerm:
        if len(options.searchTerm) > 0:
            # filter out invalid rows with no name
            df = df[pd.notnull(df['Protein names'])]

            newdf = df[ df['Protein names'].str.contains(options.searchTerm[0]) ]
            for s in options.searchTerm[1:]:
                newdf = pd.concat(
                    [
                        newdf,
                        df[ df['Protein names'].str.contains(s) ],
                    ],
                    axis = 0
                )
        else:
            newdf = df
    else:
        newdf = df

    proteins = newdf['Protein IDs'].tolist()

    # unwrap protein groups
    target_protein_ids = []
    for pg in proteins:
        pg_proteins = pg.split(';')
        target_protein_ids.extend(pg_proteins)

    save = True
    if save:
        peptides = getPeptides(target_protein_ids, options)
        getEICData(peptides, outPath, options)
#        pickle.dump(peptides, open('peptides.p','wb'))
#    else:
#        peptides = pickle.load(open('peptides.p','rb'))

    if options.profile:
        drawProfilePlots(peptides, outPath, options)

    if options.plot:
        drawPlots(peptides, outPath, options)

    return

if __name__ == '__main__':
    options =  parser.parse_args()
    main(options)
