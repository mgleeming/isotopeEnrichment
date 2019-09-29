import pandas as pd
import os, sys, pymzml, pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import scipy.optimize as opt
import shutil
import argparse
import pyteomics.mass as pymass

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
                    help = 'Text search terms used to filter proteins that are to be indcluded in the analysis. If omitted, all peptides will be included. To specify multiple search terms, include multiple argument/value pairs. For exmaple --searchTerm 60S --searchTerm 40S --searchTerm 30S'
                    )
parser.add_argument('--mzmlFile',
                    required = True,
                    type = str,
                    help = 'File path of mzML data file'
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
    def __init__(self, TIC = None, mz = None, rt = None, z = None, mods = None, mod_seq = None, sequence = None, proteinGroup = None):
        global options
        self.rts = []
        self.ints = []
        self.mzs = []
        self.mzInts = []

        self.TIC = float(TIC)
        self.mz = float(mz)
        self.rt = float(rt)
        self.z = int(z)
        self.mods = mods
        self.mod_seq = mod_seq
        self.sequence = sequence
        self.proteinGroup = proteinGroup

        self.getTargets(options)
        self.getFormula()
        return

    def getTargets(self, options):
        self.targets = []
        for isotope in range(options.clusterWidth):
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
    with open(options.msmsScansFile,'r') as if1:
        for i, l in enumerate(if1):

            # skip header row
            if i == 0: continue

            parts = l.split('\t')

            sequence = parts[11]
            peptide_proteins = parts[31]

            if len(sequence) < 2 or len(peptide_proteins) < 2 or 'CON__' in peptide_proteins: continue

            # test if this MSMS spectrum matches a protein in target_protein_ids
            match = False
            for row_protein in peptide_proteins.split(';'):
                if row_protein in target_protein_ids:
                    match = True

            if not match: continue

            # skip cases where charge in undetermined
            if int(parts[16]) == 0: continue

            p = Peptide(
                TIC = parts[4],
                mz = parts[14],
                rt = parts[2],
                z = parts[16],
                mods = parts[29],
                mod_seq = parts[30],
                sequence = parts[11],
                proteinGroup = peptide_proteins
            )

            peptides.append(p)

    peptides.sort(key = lambda x: x.TIC, reverse = True)
    return peptides

def getEICData(peptides, outPath, options):
    base = os.path.basename(options.mzmlFile)
    output = 'result_%s.dat' %base.split('.')[0]

    of1 = open(os.path.join(outPath, output),'wt')

    of1.write('Sequence, Formula, m/z, Retantion Time (min), Charge, Modifications, Protein Gorup,  %s\n'%(
        ', '.join(['Intensity %s' % str(x) for x in range(options.clusterWidth)]),
        )
    )

    # get EIC
    spectra = readSpectra(options.mzmlFile, msLevel = 1)
    for s in spectra:
        time, mzs, ints, lvl = s
        for i, p in enumerate(peptides):
            if abs(p.rt - time) < 2:
                ll = p.targets[0][0]
                hl = p.targets[0][1]

                p.rts.append(time)
                mzIndices = np.where((mzs > ll) & (mzs < hl))
                p.ints.append(np.sum(ints[mzIndices]))
                p.mzs.append(mzs)
                p.mzInts.append(ints)

    # fit gaussians
    for p in peptides:
        cen = p.rt
        wid = max(p.rts) - min(p.rts)
        amp = max(p.ints)
        p0 = [amp,cen,wid]
        try:
            RTpopt, RTpcov= opt.curve_fit(
                    lambda x, amp, cen, wid: gauss(x, amp, cen, wid, scale = 1),
                    p.rts, p.ints, p0=p0, maxfev = 2000
            )

            p.fits = [RTpopt, RTpcov]

            fit_ints = gauss(p.rts, *p.fits[0])
            index = np.argmax(fit_ints)
            p.maxrt = p.rts[index]
            p.maxrtindex = index

            abundances = []
            done = False
            for t in p.targets:
                abundances.append([])
                for index, eicRT in enumerate(p.rts):
                    if abs(index - p.maxrtindex) < options.avgNSpectra:
                        specmzs = p.mzs[index]
                        specints = p.mzInts[index]
                        mask = np.where(
                            (specmzs > t[0])
                            &
                            (specmzs < t[1])
                        )
                        intsubset = specints[mask]
                        abundances[-1].append(np.sum(intsubset))

            p.abundances = abundances
            try:
                f = p.fits
            except AttributeError:
                continue


            of1.write('%s, %s, %s, %s, %s, %s, %s, %s\n'%(
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
        except:
            pass

    of1.close()
    return

def drawPlots(peptides, outPath, opeions):
    base = os.path.basename(opeions.mzmlFile)
    output = 'result_%s.dat' % base.split('.')[0]

    figDir = 'figs_%s' % base

    path = os.path.join(outPath, figDir)
    os.makedirs(path)

    for p in peptides:
        try:
            f = p.fits
        except AttributeError:
            continue

        fig, axs = plt.subplots(2,1)
        axs[0].set_title( '%s RT: %s m/z: %s' %(p.sequence, p.rt, p.mz) )
        axs[0].set_xlabel('Retention Time (min)')
        axs[0].set_ylabel('Abundance')
        axs[0].plot(p.rts, p.ints)
        axs[0].plot(p.rts, gauss(p.rts, *p.fits[0]))

        x = list(range(len(p.abundances)))
        axs[1].bar(x, [sum(x) for x in p.abundances] , align = 'center' , alpha = 0.5)
        axs[1].set_xlabel('Isotope')
        axs[1].set_ylabel('Abundance')
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
        print(outPath)
        os.makedirs(outPath)
    except:
        print ('\nOutput directory already exists! Exiting...\n')
        sys.exit()

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

    if not options.noPlot:
        drawPlots(peptides, outPath, options)

    return

if __name__ == '__main__':
    options =  parser.parse_args()
    main(options)
