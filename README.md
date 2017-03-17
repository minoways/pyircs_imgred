# pyircs_imgred --- IRCS imaging data reduction pipeline
# Written by Yosuke Minowa (Subaru Telescope)

######### Revision History #########

####################################

*** imgred_all.py --- pipeline script ***

Usage: imgred_all.py input_list output_image [options]

Options:
  -h, --help         show this help message and exit
  --start=START      First step of the reduction sequence (default=0)
  --end=END          Last step of the reduction sequence (default=15)
  --quick            execute quick reduction (default=False)
  --flat=FLAT        Name of flat frame (default=none)
  --skyflat          Make self sky flat frame (default=False)
  --dark=DARK        Name of dark frame (default=none)
  --dist=DIST        Name of distortion database (default=none)
  --fitgeom=FITGEOM  Fitting geometry for image shift (shift|rotate,
                     default=shift)
  --fringe           Remove fringe pattern? (default=False)
  --fr_x0=FR_X0      Fringe center X coordinate (default=0.0)
  --fr_y0=FR_Y0      Fringe center Y coordinate (default=0.0)
  --fr_step=FR_STEP  Step size for fringe fitting (default=5)
  --bpm=BPM          bad pixel mask frame (default=none)
  --iternum=ITERNUM  Number of clipping iterations (default=50)
  --nsigma=NSIGMA    N-sigma reject limit (default=3.0)
  --minpix=MINPIX    minimum number of pixels in detected objects
                     (default=250)
  --hsigma=HSIGMA    sigma threshold above sky (default=3)
  --lsigma=LSIGMA    sigma threshold below sky (default=10)
  --conv=CONV        convolution kernel for objmask (default="block 3 3")
  --combine=COMBINE  type of combine operation (average|median|sum)
  --reject=REJECT    type of rejection
                     (none|minmax|ccdclip|crreject|sigclip|avsigclip|pclipor)
  --expmap=EXPMAP    exposure map file name (default=none)
  --sigmap=SIGMAP    sigma map file name (default=none)
  --erase            erase intermediate files (default=False)
  --notrace          automatic trace of object position (default=True)

Step 0: copy raw frames to the current directory
Step 1: make self sky flat, if --skyflat option is used
Step 2: flat field (and bad pixel correction, if --bpm option is used)
Step 3: distortion correction, if --dist option is used
Step 4: bright object mask
Step 5: sky subtraction
Step 6: offset measurement 
Step 7: remove unused images (bad frame etc.) 
Step 8: shift and combine
(2nd iteration for full reduction)
Step 9: bright object mask using combined image
Step 10: new object mask registration
Step 11: make self sky flat frame with new mask
Step 12: flat field with new sky flat
Step 13: distortion correction
Step 14: sky subtraction
Step 15: shift and combine
