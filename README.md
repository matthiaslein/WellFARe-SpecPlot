# WellFARe-SpecPlot
## Plot simulated UV-Vis spectra from quantum chemical data
Copyright (C) 2015  Matthias Lein

WellFAReSpecPlot: Wellington Fast Assessment of Reactions - Spectral Data Plot

positional arguments:
  file                  input file(s) with spectroscopic data
                        recognised filetypes: g09, orca

optional arguments:
  -o OUTFILE, --outfile OUTFILE
                        save plot to image file instead of displaying in gui
                        supported filetypes: .png, .ps, .eps, .pdf
                        
  --csv CSV             save to file in csv format
                        arbitrary file name (doesn't need to end in .csv)
  
  -c CUTOFF, --cutoff CUTOFF
                        cutoff value for inclusion into plots
                        Spectra contributing to less than this value (as per their Boltzmann Factor)
                        won't be shown as separate plots (default value: 0.01 [equiv 1%])
                        
  -u UPPER, --upper UPPER
                        highest frequency (in nm) for the plot
  -l LOWER, --lower LOWER
                        lowest frequency (in nm) for the plot
                        
  -b BROADENING, --broadening BROADENING
                        line broadening (in nm, default value: 3099.6)
                        
  --hwhm HWHM           half width at half peak height (only for Lorentzians;
                        in nm, default value: 7.5 nm)
                        
  --nolines             prevent printing of line spectra underneath plots
  --nonames             prevent printing of file names in plots
  --nocontr             prevent printing of contributions (Boltzmann Factors) in plots
    -t, --totalonly       only print the total plot, not individual subplots
  --flipecd             invert the handedness of the ECD data

  -f {gaussian,lorentzian}, --function {gaussian,lorentzian}
                        type of function to fit spectrum (default: gaussian)

  -p POINTS, --points POINTS
                        number of points to plot
  --colourmap {0,1,2,3}
                        choose colour map for the plot
  -v {0,1,2,3}, --verbosity {0,1,2,3}
                        increase output verbosity
