import sys
import os.path
import time

def timestamp(s):
  print (s + time.strftime("%Y/%m/%d %X"))

# ASCII FONTS from: http://patorjk.com/software/taag/
# Font = "Small"
def ProgramHeader():
  print ("#################################################################")
  print ("Wellington Fast Assessment of Reactions Spectroscopical Data Plot")
  print ("__      __   _ _ ___ _   ___     ___              ___ _     _   ")
  print ("\ \    / /__| | | __/_\ | _ \___/ __|_ __  ___ __| _ \ |___| |_ ")
  print (" \ \/\/ / -_) | | _/ _ \|   / -_)__ \ '_ \/ -_) _|  _/ / _ \  _|")
  print ("  \_/\_/\___|_|_|_/_/ \_\_|_\___|___/ .__/\___\__|_| |_\___/\__|")
  print ("                                    |_|                         ")
  print ("                                                    Version 0.01")
  print ("      WellFAReSpecPlot Copyright (C) 2015 Matthias Lein         ")
  print ("       This program comes with ABSOLUTELY NO WARRANTY           ")
  print ("        This is free software, and you are welcome to           ")
  print ("          redistribute it under certain conditions.             ")
  timestamp('Program started at: ')
  print ("#################################################################\n")

def ProgramFooter():
  print ("\n######################################################")
  print ("  ___                                ___         _    ")
  print (" | _ \_ _ ___  __ _ _ _ __ _ _ __   | __|_ _  __| |___")
  print (" |  _/ '_/ _ \/ _` | '_/ _` | '  \  | _|| ' \/ _` (_-<")
  print (" |_| |_| \___/\__, |_| \__,_|_|_|_| |___|_||_\__,_/__/")
  print ("              |___/                                   ")
  timestamp('Program terminated at: ')
  print ("#######################################################")

def ProgramAbort():
  print ("\n#################################################")
  print ("  ___                  _             _          _ ")
  print (" | _ \_  _ _ _    __ _| |__  ___ _ _| |_ ___ __| |")
  print (" |   / || | ' \  / _` | '_ \/ _ \ '_|  _/ -_) _` |")
  print (" |_|_\\\_,_|_||_| \__,_|_.__/\___/_|  \__\___\__,_|")
  timestamp('Program aborted at: ')
  print ("####################################################")
  sys.exit()
  return

def ProgramWarning(warntext=''):
  print ("\n###################################")
  print (" __      __             _           ")
  print (" \ \    / /_ _ _ _ _ _ (_)_ _  __ _ ")
  print ("  \ \/\/ / _` | '_| ' \| | ' \/ _` |")
  print ("   \_/\_/\__,_|_| |_||_|_|_||_\__, |")
  print ("                              |___/ ")
  timestamp('Warning time/date: ')
  if warntext != '':
    print ("#####################################")
    print("# ", warntext)
  print ("#####################################")
  return

def ProgramError(errortext=''):
  print ("\n#####################")
  print ("  ___                 ")
  print (" | __|_ _ _ _ ___ _ _ ")
  print (" | _|| '_| '_/ _ \ '_|")
  print (" |___|_| |_| \___/_|  ")
  timestamp('Error time/date: ')
  if errortext != '':
    print ("#######################")
    print("# ", errortext)
  print ("#######################")
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
import matplotlib.pyplot

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
    minima.append(min(i))
  return min(minima)

def findmax(listoflists):
  maxima = []
  for i in listoflists:
    maxima.append(max(i))
  return max(maxima)

# Start of the program
ProgramHeader()

# A sqrt(2) * standard deviation of 0.4 eV is 3099.6 nm. 0.1 eV is 12398.4 nm. 0.2 eV is 6199.2 nm.
stdev = 3099.6
# For Lorentzians, gamma is half bandwidth at half peak height (nm)
gamma = 12.5

if len(sys.argv) < 2:
  ProgramError("Need filename(s) as command line argument!")
  ProgramAbort()

# Create empty lists for energies, bands and osc strengths
energies = []
bands = []
strengths = []

# Excitation energies in nm
# Oscillator strengths (dimensionless)
for i in range(1,len(sys.argv)):
  if os.path.isfile(str(sys.argv[i])):
    print("Reading from file {} now.".format(str(sys.argv[i])))
    band, f, energy = extractExcitations(str(sys.argv[i]))
    bands.append(band)
    strengths.append(f)
    energies = energies + [energy]
  else:
    ProgramError("Something wrong with the file given as command line argument!")
    ProgramAbort()

# Basic check that we have the same number of bands and oscillator strengths
#if len(bands) != len(f):
#    print('Number of bands does not match the number of oscillator strengths.')
#    sys.exit()

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

print("We have spectral data from {} calculations".format(len(bands)))
for i in range(1,len(bands)+1):
  print("Calculation no {}".format(i))
  print("Gibbs energy: {}".format(energies[i-1]))
  print("Peak positions at: {}".format(bands[i-1]))
  print("Peak intensities : {}".format(strengths[i-1]))
  print("")

# Now that we know the bands, setup plot
start=np.trunc(max(findmin(bands)-50.0,0.0))
finish=np.trunc(findmax(bands)+50.0)
points=int((finish-start)*2.5)

print("Plot boundaries: {} nm and {} nm ({} points)".format(start, finish, points))

x = np.linspace(start,finish,points)

composite = 0
for i in range(0,len(bands)):
  for count,peak in enumerate(bands[i]):
      thispeak = gaussBand(x, peak, strengths[i][count], stdev)
#      thispeak = lorentzBand(x, peak, f[count], stdev, gamma)
      composite += thispeak

fig, ax = matplotlib.pyplot.subplots()
ax.plot(x,composite)
matplotlib.pyplot.xlabel('$\lambda$ / nm')
matplotlib.pyplot.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')

matplotlib.pyplot.show()

ProgramFooter()
