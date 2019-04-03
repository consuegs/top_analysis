## Requirements ##
- ROOT6
- C++11 compiler
- CMake >= 2.8 (lx-cluster: cmake28 executable)
- Boost

## Build ##
    mkdir build; cd build
    cmake ..
    make

On lx-cluster:

- use `CMSSW_7_4_X`
- use `cmake28` instead of `cmake`

## Run ##
Execute `run.x` (created in the build directory).
Without any arguments, a test module is run.
Meaningful modules can be chosen as arguments:

    $ run.x --help
    Usage: ./run.x module1[ module2[...]] [options]
    Allowed options:
      -h [ --help ]              produce help message
      -f [ --fraction ] arg (=1) Fraction of events to process
      --release                  Release mode (don't draw version labels)


## Output ##
Plots are stored as canvases in a root file.

## Adding new modules ##
To add a new module, create a new `<module>.cpp` file in `src/modules`.
The entry function has to be `extern "C" void run() {...}`.
Add `<module>` to the list in `src/modules/CMakeLists.txt`

## Changes/Additions made by the Combination ##
All scripts (in `src/modules`) used by the initial analysis are kept.
Additional scripts are denoted with a prepending `danilo_`.
The additional/changed scripts are used for the following tasks:
- `danilo_GGM_scan` : Study of kinematic variables for both GGM scans
- `danilo_SFapplication` : Based to `SFapplication`
- `danilo_acceptanceHist` : Plot signal acceptances
- `danilo_check_overlap` : Study single events, falsely removed by vetoes
- `danilo_check_scan` : Check gen particles of first GGM simulation attempts
- `danilo_cutflow_GGM` : Simple cutflows for GGM scenarios
- `danilo_datacards` : Used to create datacards
- `danilo_distributions` : Based on `distributions` with addition vetoes
- `danilo_interpolateAcc` : GGM acceptance interpolated
- `danilo_maxSensitivity` : Creates plots of most sensitive search
- `danilo_nDecays` : Gen studies of decays in GGM scenarios
- `danilo_plot1d_limit` : Plotting tool for 1D limits
- `danilo_plot2d_scan` : Plotting tool for GGM xsecs
- `danilo_plotGGMparameter` : Plotting tool for sparticle masses in GGM
- `danilo_plotOverlap` : Plotting tool for overlap between searches
- `danilo_plotSig` : Plotting tool for significance in GGM scenarios
- `danilo_plot_GGM` : Plotting tool for GGM kinematics
- `danilo_plot_postfitBKG` : Plotting tool for postfit studies
- `danilo_rValuediff` : Compare r values for different strong scenarios
- `danilo_signal_scan` : Based on `signal_scan` with additional vetoes
- `danilo_signal_scan_newTrig` : Studies for possible new trigger
- `danilo_totsignalyield` : GGM signal yields and unc.
...
