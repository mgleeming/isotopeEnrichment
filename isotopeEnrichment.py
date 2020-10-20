import os, sys, pymzml, pickle, shutil, argparse
import numpy as np
import pandas as pd
import scipy.signal
import scipy.optimize as opt
import pyteomics.mass as pymass
import matplotlib.pyplot as plt

DEFAULT_ISOTOPE_COUNTER = 6
DEFAULT_EXTRACTION_WIDTH = 0.01
DEFAULT_ISOTOPE_WEIGHT = 1
DEFAULT_SPECTRA_TO_AVERAGE = 3
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
                    type = int,
                    help = 'width (in m/z) used to produce EIC plots'
                    )
parser.add_argument('--isotopeWeight',
                    default = DEFAULT_ISOTOPE_WEIGHT,
                    type = int,
                    help = 'mass increment of isotope of interest'
                    )
parser.add_argument('--noPlot',
                    action = 'store_true',
                    help = "Don't draw peptide EIC and isotope intensity graphics"
                    )
parser.add_argument('--outDirName',
                    help = 'Results directory name',
                    default = DEFAULT_OUT_DIR_NAME
                    )

class Peptide(object):
    def __init__(self, row):
        global options

        self.fits = {k:[] for k in options.mzmlFile}
        self.rts = {k:[] for k in options.mzmlFile}
        self.ints = {k:[] for k in options.mzmlFile}
        self.mzs = {k:[] for k in options.mzmlFile}
        self.mzInts = {k:[] for k in options.mzmlFile}

        self.TIC = float(row['Total ion current'])
        self.mz = float(row['m/z'])
        self.rt = float(row['Retention time'])
        self.z = int(row['Charge'])
        self.mods = row['Modifications']
        self.mod_seq = row['Modified sequence']
        self.sequence = row['Sequence']
        self.proteinGroup = row['Proteins']

        self.getTargets(options)
        self.getFormula()

        self.abundances = {}
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

        return

    def getFormula(self):
        composition = pymass.Composition(self.sequence)
        elements = ['C', 'H', 'N', 'O', 'P', 'S']
        self.formula = ''
        for c in composition:
            self.formula += '%s%s ' %(c, composition[c])
        return

def readSpectra (mzml_file, msLevel = None):
    msrun = pymzml.run.Reader(str(mzml_file))
    for spectrum in msrun:
        lvl = spectrum['ms level']

        if msLevel:
            if lvl != msLevel: continue

        try:
            time = spectrum['scan time']
        except:
            try:
                time = spectrum['scan start time']
            except:
                print ('skipping spectrum')
                continue
        try:
            mzs = np.array(spectrum.mz, dtype = 'float64')
            ints = np.array(spectrum.i, dtype = 'uint64')
            yield time, mzs, ints, lvl
        except:
            print ('skipping spectrum')
            continue

def gauss(x, amp, cent, wid, scale = 1):
    return(amp/ (np.sqrt(2*np.pi*(wid/scale)**2 )) * np.exp(-(x-(cent/scale))**2 / (2*(wid/scale)**2)))

def getPeptides(target_protein_ids, options):
    peptides = []

    df = pd.read_csv(options.msmsScansFile, delimiter='\t')

    for index, row in df.iterrows():

        # skip unassigned or contaminant peptides
        if len(row['Sequence']) < 2 or len(row['Proteins']) < 2 or 'CON__' in row['Proteins']: continue

        # test if this MSMS spectrum matches a protein in target_protein_ids
        matches = [_ for _ in target_protein_ids if _ in row['Proteins']]
        if len(matches) == 0: continue

        # skip cases where charge in undetermined
        if int(row['Charge']) == 0: continue

        p = Peptide( row )

        peptides.append(p)

    peptides.sort(key = lambda x: x.TIC, reverse = True)
    return peptides

def getEICData(peptides, outPath, options):
#    base = os.path.basename(options.mzmlFile)
#    output = 'result_%s.dat' %base.split('.')[0]
    output = 'result.dat'

    of1 = open(os.path.join(outPath, output),'wt')

    of1.write('Sequence, Formula, m/z, Retantion Time (min), Charge, Modifications, Protein Gorup,  %s\n'%(
        ', '.join(['Intensity %s' % str(x) for x in range(options.clusterWidth)]),
        )
    )

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
        spectra = readSpectra(mzmlFile, msLevel = 1)
        for s in spectra:
            time, mzs, ints, lvl = s
            for i, p in enumerate(peptides):
                if abs(p.rt - time) < 2:
                    ll = p.targets[0][0]
                    hl = p.targets[0][1]

                    # used for picking EIC peak max
                    p.rts[mzmlFile].append(time)

                    summedEIC = 0
                    for t in p.targets:
                        mzIndices = np.where((mzs > t[0]) & (mzs < t[1]))
                        summedEIC += np.sum(ints[mzIndices])
                    p.ints[mzmlFile].append(summedEIC)

                    # used later for actual isotope abundance calculation
                    p.mzs[mzmlFile].append(mzs)
                    p.mzInts[mzmlFile].append(ints)

        # fit gaussians
        for p in peptides:
            cen = p.rt
            wid = max(p.rts[mzmlFile]) - min(p.rts[mzmlFile])
            amp = max(p.ints[mzmlFile])
            p0 = [amp,cen,wid]

            try:
                RTpopt, RTpcov= opt.curve_fit(
                        lambda x, amp, cen, wid: gauss(x, amp, cen, wid, scale = 1),
                        p.rts[mzmlFile], p.ints[mzmlFile], p0=p0, maxfev = 2000
                )
            except RuntimeError:
                # optimal fit parameters not found
                continue

            p.fits[mzmlFile] = [RTpopt, RTpcov]

            fit_ints = gauss(p.rts[mzmlFile], *p.fits[mzmlFile][0])
            index = np.argmax(fit_ints)
            p.maxrt = p.rts[mzmlFile][index]
            p.maxrtindex = index

            abundances = []
            for t in p.targets:
                abundances.append([])
                for index, eicRT in enumerate(p.rts[mzmlFile]):
                    if abs(index - p.maxrtindex) < options.avgNSpectra:
                        specmzs = p.mzs[mzmlFile][index]
                        specints = p.mzInts[mzmlFile][index]
                        mask = np.where(
                            (specmzs > t[0])
                            &
                            (specmzs < t[1])
                        )
                        intsubset = specints[mask]
                        abundances[-1].append(np.sum(intsubset))

            # used only for plotting
            p.abundances[mzmlFile] = abundances

            try:
                f = p.fits
            except AttributeError:
                continue

            of1.write('%s, %s, %s, %s, %s, %s, %s, %s, %s\n'%(
                mzmlFile,
                p.sequence,
                p.formula,
                p.mz,
                p.rt,
                p.z,
                p.mods,
                p.proteinGroup,
                ', '.join([str(sum(x)) for x in abundances]),
                )
            )

    of1.close()
    return

def drawPlots(peptides, outPath, options):

    figDir = 'figs'

    path = os.path.join(outPath, figDir)
    try:
        os.makedirs(path)
    except:
        pass

    counter = 0
    for p in peptides:

        counter += 1
        if counter > 10: break
        fig, axs = plt.subplots(len(options.mzmlFile), 2)
        for axi, mzmlFile in enumerate(options.mzmlFile):

            try:
                f = p.fits[mzmlFile][0]
            except IndexError:
                continue

            axs[axi,0].set_title( mzmlFile )
            axs[axi,0].set_xlabel('Retention Time (min)')
            axs[axi,0].set_ylabel('Abundance')
            axs[axi,0].plot(p.rts[mzmlFile], p.ints[mzmlFile])
            axs[axi,0].plot(p.rts[mzmlFile], gauss(p.rts[mzmlFile], *p.fits[mzmlFile][0]))
            axs[axi,0].set_yticklabels([])

            x = list(range(len(p.abundances[mzmlFile])))
            axs[axi,1].set_title( mzmlFile )
            axs[axi,1].bar(x, [sum(x) for x in p.abundances[mzmlFile]] , align = 'center' , alpha = 0.5)
            axs[axi,1].set_xlabel('Isotope')
            axs[axi,1].set_ylabel('Abundance')
            axs[axi,1].set_yticklabels([])

        fig.tight_layout()
        plt.savefig( os.path.join(path, '%s_rt_%s_mz_%s.png' %(p.sequence, p.rt, p.mz)))
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

    TEST = False

    if TEST:

        # get list of proteins containing '60S' string
        df = pd.read_csv(options.proteinGroupsFile, delimiter='\t')

        # filter out invalid rows with no name
        df = df[pd.notnull(df['Protein names'])]

        if len(options.searchTerm) > 0:
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

        proteins = newdf['Protein IDs'].tolist()

        # unwrap protein groups
        target_protein_ids = []
        for pg in proteins:
            pg_proteins = pg.split(';')
            target_protein_ids.extend(pg_proteins)

        peptides = getPeptides(target_protein_ids, options)

        getEICData(peptides, outPath, options)

        pickle.dump(peptides, open('save.p','wb'))
    else:
        peptides = pickle.load(open('save.p','rb'))

    if not options.noPlot:
        drawPlots(peptides, outPath, options)

    return

if __name__ == '__main__':
    options =  parser.parse_args()
    main(options)
