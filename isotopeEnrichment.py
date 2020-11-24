import os, sys, pymzml, pickle, shutil, argparse
import numpy as np
import pandas as pd
import pyteomics.mass as pymass
import matplotlib.pyplot as plt

DEFAULT_MIN_CLUSTER_WIDTH = 5
DEFAULT_INTEGRATION_WIDTH = 5
DEFAULT_EXTRACTION_LENGTH = 0.5
DEFAULT_ISOTOPE_WEIGHT = 1.00727647
DEFAULT_SPECTRA_TO_AVERAGE = 3

parser = argparse.ArgumentParser(
    description = 'Extract peptide isotope abundances from LCMS data'
)

parser.add_argument('--evidenceFile',
                    required = True,
                    type = str,
                    help = 'MaxQuant evidence.txt file'
                    )
parser.add_argument('--mzmlFileDir',
                    required = True,
                    type = str,
                    help = 'Directory containing mzML files. Files in this directory with an .mzML extention (case insensitive) will be matched.'
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
parser.add_argument('--intWidth',
                    default = DEFAULT_INTEGRATION_WIDTH,
                    type = float,
                    help = 'width (in ppm) used to integrate peaks'
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
parser.add_argument('--plot',
                    action = 'store_true',
                    help = 'Draw peptide EIC and isotope intensity graphics. This is very time consuming if the number of target peptides is large.'
                    )
parser.add_argument('--plotTopN',
                    type = int,
                    help = 'Draw peptide EIC and isotope intensity graphics for only the top N most abundant target peptides (measured by MS1 precursor intensity). Only active if the --plot flag is given'
                    )
parser.add_argument('--minSamples',
                    type = int,
                    help = 'Discard peptides observed in fewer than this number of samples',
                    )
parser.add_argument('--outDir',
                    help = 'Results directory name',
                    default = os.path.join(os.getcwd(), 'results')
                    )

class Isotope(object):

    def __init__(self, index, mzLL, mzHL):

        self.index = index
        self.mzLL = mzLL
        self.mzHL = mzHL
        self.intensityArray = []
        return

    def integrate(self, indicies):
        self.intensity = sum(self.intensityArray[indicies[0]:indicies[1]+1])
        return

class Target(object):

    def __init__(self, row, attrs, options):

        # inherit attributes from peptide instance
        for attr, val in attrs.items():

            # some are not needed though
            if attr in ['targets', 'dataFiles', 'intensities']: continue
            setattr(self, attr, val)

        self.targetID = row['id']
        self.dataFile = row['Raw file']

        # this can differ a bit between data files
        self.rt = float(row['Retention time'])
        self.intensity = int(row['Intensity']) if not np.isnan(row['Intensity']) else 0

        # used in every case
        self.rts = []
        self.isotopes = []

        # use if --plot is given
        self.mzs = []
        self.mzInts = []
        self.eicInts = []

        for isotope in range(0, self.clusterWidth):

            protonatedMass = self.mass + float(self.z) * float(options.isotopeWeight)

            isotopemass = protonatedMass + float(isotope) * float(options.isotopeWeight)

            isotopemz = isotopemass / float(self.z)

            mzLL = isotopemz - options.intWidth / 1000000 * isotopemz
            mzHL = isotopemz + options.intWidth / 1000000 * isotopemz

            self.isotopes.append( Isotope ( isotope, mzLL, mzHL ))

        self.minTargetLL = min([_.mzLL for _ in self.isotopes]) - 1
        self.maxTargetHL = max([_.mzHL for _ in self.isotopes]) + 1
        return

    def getNSpec(self, options):
        # find which spectra should be averaged
        # use MQ rt and N spectra either side

        # find closest spectrum to MQ rt
        rtdiffs = np.abs(self.rt - np.asarray(self.rts))
        self.nearest = np.argmin(rtdiffs)

        # get N spec either side of this
        self.specIndicies = [self.nearest - options.avgNSpectra, self.nearest + options.avgNSpectra]

        return

    def integrateIsotopes(self):
        for iso in self.isotopes:
            iso.integrate(self.specIndicies)
        return

    def getMonoEIC(self):
        iso = self.isotopes[0]
        for i, mzs in enumerate(self.mzs):
            ints = self.mzInts[i]
            mask = np.where( (mzs > iso.mzLL) & (mzs < iso.mzHL))
            intsubset = ints[mask]
            self.eicInts.append(np.sum(intsubset))
        return

    def addSpec(self, time, mzs, ints, options):

        self.rts.append(time)

        # integrate spectrum for each isotope and add
        for iso in self.isotopes:
            mask = np.where( (mzs > iso.mzLL) & (mzs < iso.mzHL))
            iso.intensityArray.append( np.sum(ints[mask]) )

        if options.plot:
            mask = np.where( (mzs > self.minTargetLL) & (mzs < self.maxTargetHL) )
            self.mzs.append(mzs[mask])
            self.mzInts.append(ints[mask])


class Peptide(object):

    def __init__(self, row, options):

        self.sequence = row['Sequence'].upper()
        self.modSequence = row['Modified sequence'].upper()
        self.mods = row['Modifications']

        # list of target peptide intensities for ordering peptide instances
        try:
            self.intensities = [int(row['Intensity'])]
        except ValueError: # nan val
            self.intensities = []

        # seems to be the theoretical mass of the peptide + mod
        # m/z field seems to be measured m/z
        # -- use this plus charge state to calculate target m/z
        self.mass = float(row['Mass'])
        self.z = int(row['Charge'])
        totalmass = self.mass + self.z * options.isotopeWeight
        self.mz = totalmass / self.z

        self.targets, self.dataFiles = [], []

        if options.specialResidue:
            self.specialResidueCount = sum(
                [self.sequence.count(sr.upper()) for sr in options.specialResidue]
            )
        else:
            self.specialResidueCount = 0

        # the number of isotopes that should be inspected
        self.clusterWidth = options.minClusterWidth
        if options.addSpecialResidues:
            self.clusterWidth += self.specialResidueCount

        # determine peptide formula - NB this excludes PTM contributions
        composition = pymass.Composition(self.sequence)
        elements = ['C', 'H', 'N', 'O', 'P', 'S']
        self.formula = ''
        for c in composition:
            self.formula += '%s%s ' %(c, composition[c])

        # unique key for this peptide
        self.pepKey = str(self.specialResidueCount) + '_' + self.modSequence + '_' + str(self.z)
        return

    def addTarget(self, row):

        #check if target already exists
        if row['Raw file'] in self.dataFiles:

            # duplicate samplings can still exist
            # in these cases, take most intense peptide
            # --- is there any problem with this?
            target = [t for t in self.targets if t.dataFile == row['Raw file']]
            assert len(target) == 1
            target = target[0]

            rowInt = int(row['Intensity']) if not np.isnan(row['Intensity']) else 0
            if rowInt > target.intensity:
                target.rt = float(row['Retention time'])
                target.intensity = int(row['Intensity'])

            return

        # not a duplicate sampling of an already seen peptide
        try:
            self.intensities.append(int(row['Intensity']))
        except ValueError: # nan val
            self.intensities.append(0)

        # want to pass peptide attributes to all targets
        attrs = vars(self)

        self.targets.append( Target( row, attrs, options ) )
        self.dataFiles.append(row['Raw file'])
        return

    def processTargets(self, options):
        for t in self.targets:

            t.getNSpec(options)
            t.integrateIsotopes()

            if options.plot:
                t.getMonoEIC()

def getEICData(peptides, options):

    mzmlFileDir = options.mzmlFileDir
    mzmlFiles = [os.path.join(mzmlFileDir, f) for f in os.listdir(mzmlFileDir) if '.mzml' in f.lower()]

    print('\n\tExtracting MS data')
    for mzmlFile in mzmlFiles:

        print('\t\tprocessing %s' %mzmlFile)

        # get summed EIC intensity
        spectra = pymzml.run.Reader(mzmlFile)

        for s in spectra:

            if s.ms_level != 1: continue

            time = s.scan_time_in_minutes()
            mzs = s.mz
            ints = s.i

            for p in peptides:
                for t in p.targets:

                    if abs(t.rt - time) > options.eicLength:
                        continue

                    if str(t.dataFile) != str(os.path.basename(mzmlFile).split('.')[0]):
                        continue

                    t.addSpec(time, mzs, ints, options)

    return peptides

def getPeptides(options):

    mzmlFileDir = options.mzmlFileDir
    mzmlFiles = [os.path.join(mzmlFileDir, f) for f in os.listdir(mzmlFileDir) if '.mzml' in f.lower()]

    print('\n\tDetected MzML files are:')

    for mzml in mzmlFiles:
        print('\t\t%s'%mzml)

    print('\n\tGrouping peptides')

    # read evidence.txt file
    df = pd.read_csv(options.evidenceFile, delimiter='\t')

    # filter to only non contaminant and non decoy
    df = df[df['Reverse'] != '+']
    df = df[df['Potential contaminant'] != '+']

    peptides = {}

    for index, row in df.iterrows():

        # duplicate samplings can be present
        testPepKey = row['Modified sequence'].upper() + '_ ' + str(row['Charge'])

        found = [1 for mzf in mzmlFiles if row['Raw file'] in mzf]
        if sum(found) < 1: continue

        try:
            peptides[testPepKey].addTarget(row)
        except KeyError:
            peptides[testPepKey] = Peptide(row, options)
            peptides[testPepKey].addTarget(row)

    peptides = list(peptides.values())

    # sort peptides by MQ measured intensity
    peptides.sort(key=lambda p: sum(p.intensities), reverse=True)

    print('\t\tFound %s peptides' %len(peptides))

    # discard peptides observed in too few samples
    if options.minSamples:
        print('\t\t\tRemoving peptides observed in fewer than %s of %s samples...' %(
            options.minSamples, len(mzmlFiles)
        ))

        peptides = [p for p in peptides if len(p.dataFiles) >= options.minSamples]

        print('\t\t%s peptides after filtering' %len(peptides))

    return peptides

def writeOutputTable(peptides):

    print('\n\tWriting output tables')

    # get mzx cluster width - used to add headers
    maxClusterWidth = max( [p.clusterWidth for p in peptides] )

    of1 = open(os.path.join(options.outDir, 'intensities.tsv'),'wt')

    text = '\t'.join(
        ['ID', 'File', 'Sequence', 'Special Residue Count', 'Formula', 'Mass', 'Charge', 'm/z', 'Retantion Time (min)', 'Modifications'] + ['Intensity %s' % str(x) for x in range(maxClusterWidth)]
    )

    of1.write('%s\n'%(text))
    for p in peptides:
        for t in p.targets:

            text = '\t'.join(
                [str(x) for x in [
                    t.targetID, t.dataFile, t.sequence, t.specialResidueCount, t.formula, t.mass, t.z, t.mz, t.rt, t.mods,]
                ] + [
                    str(int(iso.intensity)) for iso in t.isotopes
                ]
            )
            of1.write('%s\n'%(text))

    of1.close()
    return

def drawPlots(peptides):

    print('\n\tDrawing plots')
    for peptidei, p in enumerate(peptides):
        print('\t\tplotting %s' %peptidei)

        if options.plotTopN:
            if peptidei > options.plotTopN:
                return

        path = os.path.join(options.outDir, 'figs', p.pepKey)
        try:
            os.makedirs(path)
        except:
            pass

        for t in p.targets:

            fig, ax = plt.subplots(3, 1)
            title = str(t.sequence) + ' mz: ' + str('%.4f'%t.mz) + ' z: ' + str(t.z) + ' rt: ' + str('%.2f'%t.rt)
            fig.suptitle(title, fontsize="x-large")

            # plot chromatogram
            ax[0].plot( t.rts, t.eicInts)
            ax[0].axvline(x = t.rt, color = 'red', ls = 'dashed', lw = 0.73)
            ax[0].axvline(x = t.rts[min(t.specIndicies)], color = 'blue', ls = 'dashed', lw = 0.73)
            ax[0].axvline(x = t.rts[max(t.specIndicies)], color = 'blue', ls = 'dashed', lw = 0.73)

            # plot spectrum
            ax[1].plot( t.mzs[t.nearest], t.mzInts[t.nearest],)
            ax[1].set_xlim(( t.minTargetLL, t.maxTargetHL))

            # plot bar chart
            yvals = [i.intensity for i in t.isotopes]
            xvals = [i.index for i in t.isotopes]
            barlist = ax[2].bar( xvals, yvals , align = 'center' , alpha = 0.5 )

            # add colorings to spectral regions and bar chart
            color = plt.cm.hsv(np.linspace(0,1,len(t.isotopes)))
            for isoi, iso in enumerate(t.isotopes):
                ax[1].axvspan(iso.mzLL, iso.mzHL, color=color[isoi], alpha=0.2)
                barlist[isoi].set_color(color[isoi])

            figName = os.path.join(path, '%s_%s.png' %(t.dataFile, t.targetID))
            fig.tight_layout()
            plt.savefig(figName)
            plt.clf()
            plt.cla()
            plt.close()
    return

def main(options):

    print('Started isotopeEnrichment')

    try:
        os.makedirs(options.outDir)
    except:
        print ('\nOutput directory already exists! Exiting...\n')
        sys.exit()

    peptides = getPeptides(options)

    peptides = getEICData(peptides, options)

    for p in peptides:
        p.processTargets(options)

    writeOutputTable(peptides)

    if options.plot:
        drawPlots(peptides)

    print('All done!')
    return

if __name__ == '__main__':

    options =  parser.parse_args()
    main(options)
