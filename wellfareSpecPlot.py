import sys
import os.path
import time

def timestamp(s):
  print (s + time.strftime("%Y/%m/%d %X"))

# ASCII FONTS from: http://patorjk.com/software/taag/
# Font = "Small"
def ProgramHeader():
  print ("###############################################################################")
  print ("Wellington Fast Assessment of Reactions - Spectroscopical Data Plot")
  print ("  __      __   _ _ ___ _   ___     ___              ___ _     _    ")
  print ("  \ \    / /__| | | __/_\ | _ \___/ __|_ __  ___ __| _ \ |___| |_  ")
  print ("   \ \/\/ / -_) | | _/ _ \|   / -_)__ \ '_ \/ -_) _|  _/ / _ \  _| ")
  print ("    \_/\_/\___|_|_|_/_/ \_\_|_\___|___/ .__/\___\__|_| |_\___/\__| ")
  print ("                                    |_|                            ")
  print ("                                                       Version 0.01")
  print ("         WellFAReSpecPlot Copyright (C) 2015 Matthias Lein         ")
  print ("          This program comes with ABSOLUTELY NO WARRANTY           ")
  print ("           This is free software, and you are welcome to           ")
  print ("             redistribute it under certain conditions.             ")
  timestamp('Program started at: ')
  print ("###############################################################################\n")

def ProgramFooter():
  print ("\n###############################################################################")
  print ("  ___                                ___         _    ")
  print (" | _ \_ _ ___  __ _ _ _ __ _ _ __   | __|_ _  __| |___")
  print (" |  _/ '_/ _ \/ _` | '_/ _` | '  \  | _|| ' \/ _` (_-<")
  print (" |_| |_| \___/\__, |_| \__,_|_|_|_| |___|_||_\__,_/__/")
  print ("              |___/                                   ")
  timestamp('Program terminated at: ')
  print ("###############################################################################")

def ProgramAbort():
  print ("\n###############################################################################")
  print ("  ___                  _             _          _ ")
  print (" | _ \_  _ _ _    __ _| |__  ___ _ _| |_ ___ __| |")
  print (" |   / || | ' \  / _` | '_ \/ _ \ '_|  _/ -_) _` |")
  print (" |_|_\\\_,_|_||_| \__,_|_.__/\___/_|  \__\___\__,_|")
  timestamp('Program aborted at: ')
  print ("###############################################################################")
  sys.exit()
  return

def ProgramWarning(warntext=''):
  print ("\n###############################################################################")
  print (" __      __             _           ")
  print (" \ \    / /_ _ _ _ _ _ (_)_ _  __ _ ")
  print ("  \ \/\/ / _` | '_| ' \| | ' \/ _` |")
  print ("   \_/\_/\__,_|_| |_||_|_|_||_\__, |")
  print ("                              |___/ ")
  timestamp('Warning time/date: ')
  if warntext != '':
    print ("###############################################################################")
    print("# ", warntext)
  print ("###############################################################################")
  return

def ProgramError(errortext=''):
  print ("\n###############################################################################")
  print ("  ___                 ")
  print (" | __|_ _ _ _ ___ _ _ ")
  print (" | _|| '_| '_/ _ \ '_|")
  print (" |___|_| |_| \___/_|  ")
  timestamp('Error time/date: ')
  if errortext != '':
    print ("###############################################################################")
    print("# ", errortext)
  print ("###############################################################################")
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
  gibbsfree = 0.0
  f = open(filename,'r')
  program = "N/A"
  # Determine which QM program we're dealing with
  for line in f:
    if line.find("Entering Gaussian System, Link 0=g09") != -1:
      program = "g09"
      break
    elif line.find("* O   R   C   A *") != -1:
      program = "orca"
      break
  f.close()

  # Excitations READING SECTION
  excit = []
  # Read through Gaussian file, read *last* set of "Excitation energies and oscillator strengths:"
  if program == "g09":
    f = open(filename,'r')
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
      readBuffer=i.split()
      bands.append(float(readBuffer[6]))
      oscstr.append(float(readBuffer[8][2:]))
    f.close()
  # Read through Gaussian file, read *last* Gibbs free energy
  if program == "g09":
    gibbsfree = 0.0
    f = open(filename,'r')
    for line in f:
      if line.find("Sum of electronic and thermal Free Energies") != -1:
        gibbsfree = float(line.split()[7])
        #print("Gibbs free energy: {}".format(gibbsfree))
    f.close()
  # Read through ORCA file, read *last* set of cartesian coordinates
  elif program == "orca":
    f = open(filename,'r')
    for line in f:
      if line.find(" ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS") != -1:
        del excit[:]
        for i in range(0,4):
          readBuffer = f.__next__()
        while True:
          readBuffer = f.__next__()
          # if line is NOT empty: read results
          if readBuffer and readBuffer.strip():
            excit.append(readBuffer)
          else:
            break
    for i in excit:
      readBuffer=i.split()
      bands.append(float(readBuffer[2]))
      oscstr.append(float(readBuffer[3]))
    f.close()
  # Read through ORCA file, read *last* Gibbs free energy
  if program == "orca":
    f = open(filename,'r')
    for line in f:
      if line.find("Final Gibbs free enthalpy") != -1:
        gibbsfree = float(line.split()[5])
        #print("Gibbs free energy: {}".format(gibbsfree))
    f.close()
  #print("Excitations: ", bands, oscstr, gibbsfree)
  return bands, oscstr, gibbsfree

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
    bandshape = 1.3062974e8 * (strength / (1e7/stdev))  * np.exp(-(((1.0/x)-(1.0/band))/(1.0/stdev))**2)
    return bandshape

def lorentzBand(x, band, strength, stdev, gamma):
    "Produces a Lorentzian curve"
    bandshape = 1.3062974e8 * (strength / (1e7/stdev)) * ((gamma**2)/((x - band)**2 + gamma**2))
    return bandshape

# Start of the program
ProgramHeader()

if len(sys.argv) < 2:
  ProgramError("Need filename(s) as command line argument!")
  ProgramAbort()

# Create empty lists for energies, bands and osc strengths
energies = []
bands = []
strengths = []
names = []

# Excitation energies in nm
# Oscillator strengths (dimensionless)
for i in range(1,len(sys.argv)):
  if os.path.isfile(str(sys.argv[i])):
    print("Reading from file {} now.".format(str(sys.argv[i])))
    band, f, energy = extractExcitations(str(sys.argv[i]))
    names.append(str(sys.argv[i]))
    if band == [] or f == []:
      ProgramError("No spectral data found in this file!")
      ProgramAbort()
    if energy == []:
      ProgramError("No thermodynamic data (Gibbs free energy) found in this file!")
      ProgramAbort()
    if len(band) != len (f):
      ProgramError("Inconsistency with # of bands and # of osc strengths in this file!")
      ProgramAbort()
    bands.append(band)
    strengths.append(f)
    energies = energies + [energy]
  else:
    ProgramError("Something wrong with the file given as command line argument!")
    ProgramAbort()

# Convert absolute energies into relative energies
originalmin = min(energies)
for i in range(0,len(energies)):
  energies[i] = (energies[i] - originalmin)*627.5095

# Determine Boltzmann factors for all components
boltzmann = np.zeros(len(energies))
for i in range(0,len(energies)):
  boltzmann[i] = np.exp((-1.0*energies[i]/(298.15*0.0019872041)))

# Initialise a stdev array for the peak broadening
# A sqrt(2) * standard deviation of 0.4 eV is 3099.6 nm. 0.1 eV is 12398.4 nm. 0.2 eV is 6199.2 nm.
stdevs = np.full([len(energies),len(bands[0])], 3099.6)

# For Lorentzians, gamma is half bandwidth at half peak height (nm)
gammas = np.full([len(energies),len(bands[0])], 7.5)

print("We have spectral data from {} calculations".format(len(bands)))
for i in range(1,len(bands)+1):
  print("Data from file no {}: {}".format(i, names[i-1]))
  print("Relative Gibbs energy: {:.3f}".format(energies[i-1]))
  print("Boltzmann factor: {:.3f}".format(boltzmann[i-1]))
  print("Contribution: {:.1f}%".format((boltzmann[i-1]/np.sum(boltzmann))*100))
  print("Peak positions at: {}".format(bands[i-1]))
  print("Peak intensities : {}\n".format(strengths[i-1]))

# Now that we know the bands, setup plot
start=np.trunc(max(findmin(bands)-50.0,0.0))
finish=np.trunc(findmax(bands)+50.0)
points=int((finish-start)*2.5)

print("\nPlot boundaries: {} nm and {} nm ({} points)".format(start, finish, points))

x = np.linspace(start,finish,points)

individual = []
composite = 0
for i in range(0,len(bands)):
  individual.append(0)
  for count,peak in enumerate(bands[i]):
      thispeak = (boltzmann[i]/np.sum(boltzmann))*gaussBand(x, peak, strengths[i][count], stdevs[i][count])
#      thispeak = (boltzmann[i]/np.sum(boltzmann))*lorentzBand(x, peak, strengths[i][count], stdevs[i][count], gammas[i][count])
      composite += thispeak
      individual[i] += thispeak

# Find out how many structures actually contribute significantly (>1%)
sigstruct = 0
for i in range(0,len(bands)):
  if (boltzmann[i]/np.sum(boltzmann)) > 0.01:
    sigstruct += 1

print("\nThere are {} structures that contribute significantly (>1%)".format(sigstruct))

#colourmap = plt.cm.Spectral(np.linspace(0, 1, len(bands)))
#colourmap = plt.cm.rainbow(np.linspace(0, 1, len(bands)))
colourmap = plt.cm.gnuplot(np.linspace(0, 1, len(bands)))
#colourmap = plt.cm.seismic(np.linspace(0, 1, len(bands)))

# Setup for composite plot and each significantly contr. structure
fig, ax = plt.subplots(nrows=sigstruct+1,sharex=True,sharey=False)
ax[0].plot(x,composite)

# Go through all energies in order, but index in variable "count"
for count, i in enumerate(np.argsort(energies)):
  if (boltzmann[i]/np.sum(boltzmann)) > 0.01:
    ax[count+1].plot(x,individual[i],color=colourmap[i])
    ax[count+1].text(0.8, 0.5,'{}\n Contribution: {:.1f}%'.format(names[i],(boltzmann[i]/np.sum(boltzmann))*100),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax[count+1].transAxes)
  ax[0].plot(x,individual[i],color=colourmap[i],linestyle='--')
ax[0].text(0.8, 0.5,'All contributions',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax[0].transAxes)
plt.xlabel('$\lambda$ / nm')
plt.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')

plt.show()

ProgramFooter()
