import sys
import os.path
import time


def timestamp(s):
    print(s + time.strftime("%Y/%m/%d %X"))


# ASCII FONTS from: http://patorjk.com/software/taag/
# Font = "Small"
def ProgramHeader():
    print("###############################################################################")
    print("Wellington Fast Assessment of Reactions - Spectroscopical Data Plot")
    print("  __      __   _ _ ___ _   ___     ___              ___ _     _    ")
    print("  \ \    / /__| | | __/_\ | _ \___/ __|_ __  ___ __| _ \ |___| |_  ")
    print("   \ \/\/ / -_) | | _/ _ \|   / -_)__ \ '_ \/ -_) _|  _/ / _ \  _| ")
    print("    \_/\_/\___|_|_|_/_/ \_\_|_\___|___/ .__/\___\__|_| |_\___/\__| ")
    print("                                      |_|                          ")
    print("                                                        Version 0.9")
    print("         WellFAReSpecPlot Copyright (C) 2015 Matthias Lein         ")
    print("          This program comes with ABSOLUTELY NO WARRANTY           ")
    print("           This is free software, and you are welcome to           ")
    print("             redistribute it under certain conditions.             ")
    timestamp('Program started at: ')
    print("###############################################################################\n")


def ProgramFooter():
    print("\n###############################################################################")
    print("  ___                                ___         _    ")
    print(" | _ \_ _ ___  __ _ _ _ __ _ _ __   | __|_ _  __| |___")
    print(" |  _/ '_/ _ \/ _` | '_/ _` | '  \  | _|| ' \/ _` (_-<")
    print(" |_| |_| \___/\__, |_| \__,_|_|_|_| |___|_||_\__,_/__/")
    print("              |___/                                   ")
    timestamp('Program terminated at: ')
    print("###############################################################################")


def ProgramAbort():
    print("\n###############################################################################")
    print("  ___                  _             _          _ ")
    print(" | _ \_  _ _ _    __ _| |__  ___ _ _| |_ ___ __| |")
    print(" |   / || | ' \  / _` | '_ \/ _ \ '_|  _/ -_) _` |")
    print(" |_|_\\\_,_|_||_| \__,_|_.__/\___/_|  \__\___\__,_|")
    timestamp('Program aborted at: ')
    print("###############################################################################")
    sys.exit()
    return


def ProgramWarning(warntext=''):
    print("\n###############################################################################")
    print(" __      __             _           ")
    print(" \ \    / /_ _ _ _ _ _ (_)_ _  __ _ ")
    print("  \ \/\/ / _` | '_| ' \| | ' \/ _` |")
    print("   \_/\_/\__,_|_| |_||_|_|_||_\__, |")
    print("                              |___/ ")
    timestamp('Warning time/date: ')
    if warntext != '':
        print("###############################################################################")
        print("# ", warntext)
    print("###############################################################################\n")
    return


def ProgramError(errortext=''):
    print("\n###############################################################################")
    print("  ___                 ")
    print(" | __|_ _ _ _ ___ _ _ ")
    print(" | _|| '_| '_/ _ \ '_|")
    print(" |___|_| |_| \___/_|  ")
    timestamp('Error time/date: ')
    if errortext != '':
        print("###############################################################################")
        print("# ", errortext)
    print("###############################################################################")
    return


# Check for numpy and matplotlib, try to exit gracefully if not found
import imp

try:
    imp.find_module('numpy')
    foundnp = True
except ImportError:
    foundnp = False
try:
    imp.find_module('matplotlib')
    foundplot = True
except ImportError:
    foundplot = False
if not foundnp:
    ProgramError("Numpy is required. Exiting")
    ProgramAbort()
if not foundplot:
    ProgramError("Matplotlib is required. Exiting")
    ProgramAbort()
import numpy as np
import matplotlib.pyplot as plt


def extractExcitations(filename):
    bands = []
    oscstr = []
    ecdstr = []
    gibbsfree = 0.0
    try:
        f = open(filename, 'r')
    except:
        ProgramWarning("Can't open file {}".format(filename))
        return bands, oscstr, gibbsfree, ecdstr
    program = "N/A"
    # Determine which QM program we're dealing with
    for line in f:
        if line.find("Entering Gaussian System, Link 0=g09") != -1:
            if args.verbosity >= 3:
                print("{} is a Gaussian file".format(filename))
            program = "g09"
            break
        elif line.find("* O   R   C   A *") != -1:
            if args.verbosity >= 3:
                print("{} is an Orca file".format(filename))
            program = "orca"
            break
    f.close()

    # Excitations READING SECTION
    excit = []
    # Read through Gaussian file, read *last* set of "Excitation energies and oscillator strengths:"
    if program == "g09":
        f = open(filename, 'r')
        for line in f:
            if line.find("Excitation energies and oscillator strengths:") != -1:
                del excit[:]
                while True:
                    readBuffer = f.__next__()
                    if readBuffer.find("Excited State") != -1:
                        excit.append(readBuffer)
                    elif readBuffer.find("Leave Link") != -1:
                        break
        for i in excit:
            readBuffer = i.split()
            bands.append(float(readBuffer[6]))
            oscstr.append(float(readBuffer[8][2:]))
        f.close()
        # Read ECD Data from Gaussian file
        f = open(filename, 'r')
        for line in f:
            if line.find("R(velocity)    E-M Angle") != -1:
                del excit[:]
                while True:
                    readBuffer = f.__next__()
                    if readBuffer.find("del") != -1:
                        break
                    else:
                        excit.append(readBuffer)
        ecdstr = []
        for i in excit:
            readBuffer = i.split()
            ecdstr.append(float(readBuffer[4]))
        f.close()
    # Read through Gaussian file, read *last* Gibbs free energy
    if program == "g09":
        gibbsfree = 0.0
        f = open(filename, 'r')
        for line in f:
            if line.find("Sum of electronic and thermal Free Energies") != -1:
                gibbsfree = float(line.split()[7])
                # print("Gibbs free energy: {}".format(gibbsfree))
        f.close()
    # Read through ORCA file, read *last* set of cartesian coordinates
    elif program == "orca":
        f = open(filename, 'r')
        for line in f:
            if line.find(" ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS") != -1:
                del excit[:]
                for i in range(0, 4):
                    readBuffer = f.__next__()
                while True:
                    readBuffer = f.__next__()
                    # if line is NOT empty: read results
                    if readBuffer and readBuffer.strip():
                        excit.append(readBuffer)
                    else:
                        break
        for i in excit:
            readBuffer = i.split()
            bands.append(float(readBuffer[2]))
            oscstr.append(float(readBuffer[3]))
        f.close()
    # Read through ORCA file, read *last* Gibbs free energy
    if program == "orca":
        f = open(filename, 'r')
        for line in f:
            if line.find("Final Gibbs free enthalpy") != -1:
                gibbsfree = float(line.split()[5])
                # print("Gibbs free energy: {}".format(gibbsfree))
        f.close()
    # print("Excitations: ", bands, oscstr, gibbsfree)
    return bands, oscstr, gibbsfree, ecdstr


def findmin(listoflists):
    minima = []
    for i in listoflists:
        if i != []:
            minima.append(min(i))
    return min(minima)


def findmax(listoflists):
    maxima = []
    for i in listoflists:
        if i != []:
            maxima.append(max(i))
    return max(maxima)


# Information on producing spectral curves (Gaussian and Lorentzian) is adapted from:
# P. J. Stephens, N. Harada, Chirality 22, 229 (2010).
# Gaussian curves are often a better fit for UV/Vis.
def gaussBand(x, band, strength, stdev):
    "Produces a Gaussian curve"
    bandshape = 1.3062974e8 * (strength / (1e7 / stdev)) * np.exp(-(((1.0 / x) - (1.0 / band)) / (1.0 / stdev)) ** 2)
    return bandshape


def lorentzBand(x, band, strength, stdev, gamma):
    "Produces a Lorentzian curve"
    bandshape = 1.3062974e8 * (strength / (1e7 / stdev)) * ((gamma ** 2) / ((x - band) ** 2 + gamma ** 2))
    return bandshape


# Start of the program
import argparse

parser = argparse.ArgumentParser(
    description="WellFAReSpecPlot: Wellington Fast Assessment of Reactions - Spectral Data Plot",
    epilog="recognised filetypes: g09, orca")
parser.add_argument("files", metavar='file', help="input file(s) with spectroscopic data", nargs='+',
                    default="reactant.log")
parser.add_argument("-c", "--cutoff", help="cutoff value for inclusion into plots",
                    default=0.01,
                    type=float)
parser.add_argument("-u", "--upper", help="highest frequency (in nm) for the plot", type=float)
parser.add_argument("-l", "--lower", help="lowest frequency (in nm) for the plot", type=float)
parser.add_argument("-b", "--broadening", help="line broadening (in nm)", type=float, default=3099.6)
parser.add_argument("--hwhm", help="half width at half peak height (only for Lorentzians; in nm)", type=float,
                    default=7.5)
parser.add_argument("--nolines", help="prevent printing of line spectra underneath main plots", action='store_true')
parser.add_argument("--nonames", help="prevent printing of file names in plots", action='store_true')
parser.add_argument("-t", "--totalonly", help="only print the total plot, not individual subplots", action='store_true')
parser.add_argument("--flipecd", help="invert the handedness of the ECD data", action='store_true')
parser.add_argument("-f", "--function", help="type of function to fit spectrum", choices=["gaussian", "lorentzian"],
                    default="gaussian")
parser.add_argument("-p", "--points", help="number of points to plot", type=float)
parser.add_argument("--colourmap", help="choose colour map for the plot", type=int, choices=[0, 1, 2, 3], default=0)
parser.add_argument("-v", "--verbosity", help="increase output verbosity", type=int, choices=[0, 1, 2, 3], default=1)

args = parser.parse_args()

if args.verbosity >= 2:
    ProgramHeader()

# Create empty lists for energies, bands and osc strengths
energies = []
bands = []
strengths = []
ecds = []
names = []

# Excitation energies in nm
# Oscillator strengths (dimensionless)
for i in range(0, len(args.files)):
    if os.path.isfile(args.files[i]):
        if args.verbosity >= 2:
            print("Reading from file {} now.".format(args.files[i]))
        band, f, energy, ecd = extractExcitations(args.files[i])
        if band == [] or f == [] or ecd == []:
            ProgramWarning("No spectral data found in this file!")
        elif energy == []:
            ProgramWarning("No thermodynamic data (Gibbs free energy) found in this file!")
        elif len(band) != len(f) or len(band) != len(ecd):
            ProgramWarning("Inconsistency with # of bands and # of osc strengths/ecd in this file!")
        else:
            names.append(args.files[i])
            bands.append(band)
            strengths.append(f)
            ecds.append(ecd)
            energies = energies + [energy]
    else:
        ProgramWarning("The file {} doesn't exist or is not a file".format(args.files[i]))

if len(energies) > 1:
    # Convert absolute energies into relative energies
    originalmin = min(energies)
    for i in range(0, len(energies)):
        energies[i] = (energies[i] - originalmin) * 627.5095
    # Determine Boltzmann factors for all components
    boltzmann = np.zeros(len(energies))
    for i in range(0, len(energies)):
        boltzmann[i] = np.exp((-1.0 * energies[i] / (298.15 * 0.0019872041)))
else:
    if len(energies) == 1:
        boltzmann = np.ones(len(energies))
    else:
        ProgramError("No spectral data for plotting")
        ProgramAbort()


# Initialise a stdev array for the peak broadening
# A sqrt(2) * standard deviation of 0.4 eV is 3099.6 nm. 0.1 eV is 12398.4 nm. 0.2 eV is 6199.2 nm.
stdevs = np.full([len(energies), len(bands[0])], args.broadening)

# For Lorentzians, gamma is half bandwidth at half peak height (nm)
gammas = np.full([len(energies), len(bands[0])], args.hwhm)


# Find out how many structures actually contribute significantly to the UV-Vis
sigstruct = 0
if len(names) > 1:
    for i in range(0, len(bands)):
        if (boltzmann[i] / np.sum(boltzmann)) > args.cutoff:
            sigstruct += 1
else:
    sigstruct = 1

# Flip the ECD spectrum (to show the "other" enantiomer)
if args.flipecd == True:
    ecds = np.multiply(ecds, -1.0)

# Go through all significant structures again and check if they have an ECD spectrum
ecd_sigstruct = 0
for count, i in enumerate(np.argsort(energies)):
    if (boltzmann[i] / np.sum(boltzmann)) > args.cutoff and max(np.absolute(ecds[i])) > 0.0:
        ecd_sigstruct += 1

# Now that we know the bands, setup plot
if args.lower == None:
    start = np.trunc(max(findmin(bands) - 50.0, 0.0))
else:
    start = args.lower
if args.upper == None:
    finish = np.trunc(findmax(bands) + 50.0)
else:
    finish = args.upper
if args.points == None:
    points = int((finish - start) * 2.5)
else:
    points = args.points

x = np.linspace(start, finish, points)

if args.verbosity >= 2:
    print("")
    for i in range(1, len(bands) + 1):
        print("Data from file no {}: {}".format(i, names[i - 1]))
        print("Relative Gibbs energy: {:.3f}".format(energies[i - 1]))
        print("Boltzmann factor: {:.3f}".format(boltzmann[i - 1]))
        print("Contribution: {:.1f}%".format((boltzmann[i - 1] / np.sum(boltzmann)) * 100))
        if args.verbosity >= 3:
            print("  nm     UV-Vis     ECD")
            for j in range(0, len(bands[i - 1])):
                print(" {:.1f}  {:7.5f} {:-10.5f}".format(bands[i - 1][j], strengths[i - 1][j], ecds[i - 1][j]))
        print("")
if args.verbosity >= 1:
    print("Examined {} file(s)".format(len(args.files)))
    print("Found spectral data in {} file(s)".format(len(bands)))
if args.verbosity >= 2:
    if args.function == "lorentzian":
        print("Using Lorentzians with a line broadening of {:.1f} nm and a HWHM of {:.1f} nm for plotting".format(
            args.broadening, args.hwhm))
    else:
        print("Using Gaussians with a line broadening of {:.1f} nm for plotting".format(args.broadening))
    if args.totalonly == True:
        print("Only plotting overall UV-Vis plot")
    else:
        print("Plotting {} structure(s) that contribute to >{:.1f}% to the UV-Vis spectrum".format(sigstruct,
                                                                                                   args.cutoff * 100))
    if ecd_sigstruct > 0 and args.totalonly == False:
        print("Plotting {} contributing structure(s) with ECD data".format(ecd_sigstruct))
    if ecd_sigstruct > 0 and args.totalonly == True:
        print("Only plotting overall ECD plot. {} structures with significant ECD data.".format(ecd_sigstruct))
    else:
        if args.verbosity >= 3:
            print("No ECD spectra available or no significant contribution to the spectrum")
    if ecd_sigstruct > 0 and args.flipecd == True:
        print("The ECD data has been inverted (to show the other enantiomer)")

    if args.verbosity >= 3 and sigstruct > 1:
        print("Note that the overall UV-Vis and ECD spectra *always* contain *all* contributions.")
    print("Plotting data from {} nm to {} nm ({} points)".format(start, finish, points))
    print("")

# Calculate composite spectrum and individual spectra for all UV-Vis data
individual = []
composite = 0
for i in range(0, len(bands)):
    individual.append(0)
    for count, peak in enumerate(bands[i]):
        if args.function == "lorentzian":
            thispeak = (boltzmann[i] / np.sum(boltzmann)) * lorentzBand(x, peak, strengths[i][count], stdevs[i][count],
                                                                        gammas[i][count])
        else:
            thispeak = (boltzmann[i] / np.sum(boltzmann)) * gaussBand(x, peak, strengths[i][count], stdevs[i][count])
        composite += thispeak
        individual[i] += thispeak

# Calculate composite spectrum and individual spectra for all ECD data
individual_ecd = []
composite_ecd = 0
for i in range(0, len(bands)):
    individual_ecd.append(0)
    for count, peak in enumerate(bands[i]):
        if args.function == "lorentzian":
            thispeak = (boltzmann[i] / np.sum(boltzmann)) * lorentzBand(x, peak, ecds[i][count], stdevs[i][count],
                                                                        gammas[i][count])
        else:
            thispeak = (boltzmann[i] / np.sum(boltzmann)) * gaussBand(x, peak, ecds[i][count], stdevs[i][count])
        composite_ecd += thispeak
        individual_ecd[i] += thispeak

if args.colourmap == 0:
    colourmap = plt.cm.gnuplot(np.linspace(0, 1, len(bands)))
elif args.colourmap == 1:
    colourmap = plt.cm.Spectral(np.linspace(0, 1, len(bands)))
elif args.colourmap == 2:
    colourmap = plt.cm.rainbow(np.linspace(0, 1, len(bands)))
elif args.colourmap == 3:
    colourmap = plt.cm.seismic(np.linspace(0, 1, len(bands)))

if sigstruct == 1:
    # Setup for one individual plot if there is only one
    fig, ax = plt.subplots(nrows=1, sharex=True, sharey=False)
    ax.plot(x, composite)
    ax.set_title("UV-Vis")
    stretchfactor = 1 / max(strengths[0])
    if args.nolines != True:
        for j in range(0, len(bands[0])):
            ax.vlines(bands[0][j], 0.0, ax.get_ylim()[1] * stretchfactor * strengths[0][j])
    if args.nonames != True:
        ax.text(0.8, 0.8,
                '{}'.format(names[0]),
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    plt.xlabel('$\lambda$ / nm')
    plt.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')
elif args.totalonly == True:
    # Setup for one individual plot if especially requested
    fig, ax = plt.subplots(nrows=1, sharex=True, sharey=False)
    ax.plot(x, composite)
    ax.set_title("UV-Vis")
    # add individual contributing plots underneath
    for count, i in enumerate(np.argsort(energies)):
        ax.plot(x, individual[i], color=colourmap[i], linestyle='--')
    plt.xlabel('$\lambda$ / nm')
    plt.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')
else:
    # Setup for composite plot and each significantly contr. structure for UV-Vis
    fig, ax = plt.subplots(nrows=sigstruct + 1, sharex=True, sharey=False)
    ax[0].plot(x, composite)

    # Go through all energies in order, but index in variable "count" for UV-Vis
    ax[0].set_title("UV-Vis")
    for count, i in enumerate(np.argsort(energies)):
        if (boltzmann[i] / np.sum(boltzmann)) > args.cutoff:
            ax[count + 1].plot(x, individual[i], color=colourmap[i])
            # Ensure that the y-axis starts at zero
            ax[count + 1].axis(ymin=0.0)
            if args.nonames != True:
                ax[count + 1].text(0.8, 0.5,
                                   '{}\n Contribution: {:.1f}%'.format(names[i],
                                                                       (boltzmann[i] / np.sum(boltzmann)) * 100),
                                   horizontalalignment='center', verticalalignment='center',
                                   transform=ax[count + 1].transAxes)
            else:
                ax[count + 1].text(0.8, 0.5,
                                   'Contribution: {:.1f}%'.format((boltzmann[i] / np.sum(boltzmann)) * 100),
                                   horizontalalignment='center', verticalalignment='center',
                                   transform=ax[count + 1].transAxes)
            # print("Strongest transition in structure {}: {}".format(i,max(strengths[i])))
            stretchfactor = 1 / max(strengths[i])
            if args.nolines != True:
                for j in range(0, len(bands[0])):
                    # Print vertical line spectrum scaled to 90% of the size of the y-axis (ax[count+1].get_ylim()[1])
                    ax[count + 1].vlines(bands[i][j], 0.0,
                                         0.9 * ax[count + 1].get_ylim()[1] * stretchfactor * strengths[i][j])
                    # print("Plotting band of molecule {} at: {}".format(i,bands[i][j]))
        ax[0].plot(x, individual[i], color=colourmap[i], linestyle='--')
    ax[0].text(0.8, 0.5, 'All contributions', horizontalalignment='center', verticalalignment='center',
               transform=ax[0].transAxes)
    plt.xlabel('$\lambda$ / nm')
    plt.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')

# Setup for composite plot and each significantly contr. structure for ECD
if ecd_sigstruct == 1:
    # find out which structure it is that is contributing
    for count, i in enumerate(np.argsort(energies)):
        if (boltzmann[i] / np.sum(boltzmann)) > args.cutoff and max(np.absolute(ecds[i])) > 0.0:
            ecd_struct = i
    fig, ay = plt.subplots(nrows=1, sharex=True, sharey=False)
    ay.plot(x, composite_ecd)
    ay.axhline()
    if sigstruct > 1:
        if args.nonames != True:
            ay.text(0.8, 0.8,
                    '{}\n Contribution: {:.1f}%'.format(names[ecd_struct],
                                                        (boltzmann[ecd_struct] / np.sum(boltzmann)) * 100),
                    horizontalalignment='center', verticalalignment='center', transform=ay.transAxes)
        else:
            ay.text(0.8, 0.8,
                    'Contribution: {:.1f}%'.format(
                        (boltzmann[ecd_struct] / np.sum(boltzmann)) * 100),
                    horizontalalignment='center', verticalalignment='center', transform=ay.transAxes)
    else:
        if args.nonames != True:
            ay.text(0.8, 0.8,
                    '{}'.format(names[ecd_struct]),
                    horizontalalignment='center', verticalalignment='center', transform=ay.transAxes)
    stretchfactor = 1 / max(ecds[ecd_struct])
    if args.nolines != True:
        for j in range(0, len(bands[ecd_struct])):
            ay.vlines(bands[ecd_struct][j], 0.0, ay.get_ylim()[1] * stretchfactor * ecds[ecd_struct][j])
    ay.set_title("ECD")
    plt.xlabel('$\lambda$ / nm')
    plt.ylabel('intensity / arbitrary units')
elif args.totalonly == True:
    # setup plot
    fig, ay = plt.subplots(nrows=1, sharex=True, sharey=False)
    ay.plot(x, composite_ecd)
    ay.axhline()
    # find out which structure it is that is contributing
    for count, i in enumerate(np.argsort(energies)):
        if (boltzmann[i] / np.sum(boltzmann)) > args.cutoff and max(np.absolute(ecds[i])) > 0.0:
            ay.plot(x, individual_ecd[i], color=colourmap[i], linestyle='--')
    ay.set_title("ECD")
    plt.xlabel('$\lambda$ / nm')
    plt.ylabel('intensity / arbitrary units')
elif ecd_sigstruct > 1:
    # Go through all energies in order, but index in variable "count" for ECD
    fig, ay = plt.subplots(nrows=ecd_sigstruct + 1, sharex=True, sharey=False)
    ay[0].plot(x, composite_ecd)
    ay[0].axhline()
    countpanels = 1
    for count, i in enumerate(np.argsort(energies)):
        if (boltzmann[i] / np.sum(boltzmann)) > args.cutoff and max(np.absolute(ecds[i])) > 0.0:
            ay[countpanels].plot(x, individual_ecd[i], color=colourmap[i])
            ay[countpanels].axhline()
            if args.nonames != True:
                ay[countpanels].text(0.8, 0.8, '{}\n Contribution: {:.1f}%'.format(names[i], (
                    boltzmann[i] / np.sum(boltzmann)) * 100), horizontalalignment='center', verticalalignment='center',
                                     transform=ay[countpanels].transAxes)
            else:
                ay[countpanels].text(0.8, 0.8, 'Contribution: {:.1f}%'.format((
                                                                                  boltzmann[i] / np.sum(
                                                                                      boltzmann)) * 100),
                                     horizontalalignment='center', verticalalignment='center',
                                     transform=ay[countpanels].transAxes)

            if args.nolines != True:
                stretchfactor = 1 / max(ecds[i])
                for j in range(0, len(bands[0])):
                    # Print vertical line spectrum scaled to 90% of the size of the y-axis (ax[count+1].get_ylim()[1])
                    if ecds[i][j] > 0.0:
                        ay[countpanels].vlines(bands[i][j], 0.0,
                                               0.9 * ay[countpanels].get_ylim()[1] * stretchfactor * ecds[i][j])
                    if ecds[i][j] < 0.0:
                        ay[countpanels].vlines(bands[i][j], 0.0,
                                               -0.9 * ay[countpanels].get_ylim()[0] * stretchfactor * ecds[i][
                                                   j])
                        # print("Plotting band of molecule {} at: {}".format(i,bands[i][j]))
            countpanels += 1
        ay[0].plot(x, individual_ecd[i], color=colourmap[i], linestyle='--')
    ay[0].text(0.8, 0.8, 'All contributions', horizontalalignment='center', verticalalignment='center',
               transform=ay[0].transAxes)
    ay[0].set_title("ECD")
    plt.xlabel('$\lambda$ / nm')
    plt.ylabel('intensity / arbitrary units')

plt.show()

if args.verbosity >= 2:
    ProgramFooter()
