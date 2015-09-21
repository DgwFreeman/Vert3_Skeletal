This file is written in [Markdown](http://daringfireball.net/projects/markdown/), a simple presentation format that is rendered to HTML. The format is intended to be human-readable text even without rendering, and can be opened fine in Notepad, MATLAB's editor, whatever.

# Vert3_Simple

Vert3_Simple is a simulation of muscle contraction on the sub-sarcomere level with emphasis on:

* Calcium regulatory units: the behavior of regulatory units (RUs) on the thin filament is explicitly represented
* Cooperativity: interactions between RUs and cross-bridges (XBs) are a close part of the model
* Spatially explicit modeling: XBs and RUs are given Cartesian positions that affect their force production

More information on the nonlinear kinetics of the model can be found in Tanner et al. (2012), "Filament Compliance Influences Cooperative Activation of Thin Filaments and the Dynamics of Force Production in Skeletal Muscle", PLoS Comp. Biol. 8.

## Operation

### Parallel runs

The system is set up to allow doing multiple runs of the same parameters simultaneously on multiple CPU cores. This does not need to be specially enabled in the code; it simply uses 'parfor' (see MATLAB documentation for more information) and will run in parallel as long as MATLAB is aware of multiple CPU cores. To make MATLAB aware of multiple cores, use "matlabpool('open')" (or "parpool" in more recent MATLAB versions). The number of cores MATLAB is currently aware of can be found with "matlabpool('size')".

### Main entry points

#### Multiple

Multiple is a simple GUI for running simulations. It can be started just by entering the command "multiple" at the MATLAB prompt. More information on starting GUIDE GUIs and other methods of doing so can be found in MATLAB's help, under "Creating Graphical User Interfaces / Creating GUIs with GUIDE / Saving and Running a GUIDE GUI / Running a GUI".

Multiple only lets you specify the most common parameters: RU density, XB density, pCa, the XB stiffness scaling constant, the sarcomere length, and the titin filament length. (See below for more information on parameters.) Values can be selected from the listboxes. Multiple values can be selected, as multiple pCa levels are to begin with, and this results in all possible combinations being run. For example, if random RU densities of 1.0 and 0.7, and uniform XB densities of 0.6 and 0.5 are selected, four different simulations will be run, with RU=1.0 and XB=0.6, RU=1.0 and XB=0.5, etc.

The # of Runs box controls how many simulation runs are performed for each combination of parameters. Runs are done in parallel if possible (see below). If # of Runs is set to zero (the default) an appropriate number of runs will be automatically calculated via SimRuns.

Additional values for RU and XB densities can be added by typing a value in the relevant box, clicking "add", and selecting the newly added values.

Once you have chosen the parameter sets that should be simulated, just click "Go". As of this writing there is no progress indicator, but you can watch folders be created.

If you want to change other parameters not covered by the GUI, it is recommended you use an alternate entry point. Failing that, you can edit the init_params file to change most other parameters of the simulation.

Data files are stored in the current directory's subdirectory 'DataFiles'. Note that if DataFiles does not already exist an error will occur. This directory is hardcoded for the time being, unfortunately. The new data will be stored in a new folder DataFiles/MultipleNNN, where NNN is a number that doesn't previously exist. For instance, if Multiple1, Multiple2, and Multiple4 all exist, the new folder will be Multiple5. See directions notes on data storage below for possible replacements of this mechanism.

#### RunSeveral

RunSeveral is the main programmatic interface to running simulations.

Example usage:

    init_params;
    [Steps, Means, Vars, IndexThalf, Binder] = RunSeveral(RequestedRuns, DataParams, StartLength, pCa, StiffScale, filaments, knockout, coop, TFRateScale, tcparam);

All parameters to RunSeveral are structs, except for RequestedRuns, StartLength, and pCa. The nature of these structs is explained in the section on init_params.

RunSeveral runs several simulations under identical parameters, including the pCa. To run multiple parameter sets, you'll need to run RunSeveral several times, as by using the GUI.

As with Multiple, if RequestedRuns is zero an appropriate number will be automatically calculated via SimRuns.

RunSeveral does not handle data storage (i.e. writing to disk). Instead, it returns all data in five output matrices. These matrices were essentially chosen just to shrink the number of output parameters and should not be taken as necessary. They are defined as follows:

    Binder = Big_Binder_3;
    Steps = [MFvec; AFvec; FractXB1; FractXB2; FractCa0; FractCa1; FractCa2; ATPuse];
    Means = [MFMean; AFMean; FractXB1Mean; FractXB2Mean; FractCa0Mean; FractCa1Mean; FractCa2Mean; ATPuseMean; FractBoundMean];
    Vars = [MFVar; AFVar; FractXB1Var; FractXB2Var; FractCa0Var; FractCa1Var; FractCa2Var; ATPuseVar; FractBoundVar];
    IndexThalf = [IndexThalfMFvec; IndexThalfAFvec; IndexThalfFractXB1; IndexThalfFractXB2; IndexThalfFractCa0; IndexThalfFractCa1; IndexThalfFractCa2; IndexThalfATPuse; IndexThalfFractBound];

#### WriteText

This is a simple data storage function. It takes as arguments the appropriately named output values from RunSeveral, as well as an output directory to write files to, a pCa (used for filenames), and the timestep the simulation ran in. Note that the output directory needs to end with a file separator to indicate that it is indeed a directory.

Example usage:

    WriteText(['DataFiles' filesep 'Run683' filesep], 4, .001, binder, Steps, Means, Vars, IndexThalf);

This function mimicks the previously used output format. All files are tabular text files suitable for Excel import or Matlab import (with importdata). The files output are as follows, where #.## indicates the pCa:

*  TimeSeriesAvg\_pCa\_#.##.txt: Includes data from every time step. Values recorded are the forces produced by both filaments, ATP hydrolyzed in that timestep, the total fraction of bound XBs, and the fractions of XBs and RUs in each state. These values are averages across all runs - e.g., if two runs were carried out, and one hydrolyzed 6 ATP four milliseconds in and the other hydrolyzed 5 ATP four milliseconds in, 5.5 would be recorded.
*  SSData\_pCa\_#.##.txt: Data from the steady state of muscle contraction (see init_params for definition). The same values as in TimeSeriesAvg are recorded, but per-run, and in the form of means and variances over the steady state period.
*  HalfTimeData\_pCa\_#.##.txt: Records at which time step each value (as above) reached half of its steady-state value.
*  XBBindingData\_pCa\_#.##.txt: Records several dozen variables related to XB activity. See file headings for particulars. This is easily the largest file written out, with thousands of data points.

#### init_params

This script file loads default parameters for simulations. To run simulations under other parameters, the current way is to run init_params anyway, and then change which parameters are different.

The simulation takes many dozens of parameters. For convenience, most of these are broken up into a few structs. These structs should not be considered set in stone.

Parameters are as follows:

*  StartLength: Length of the half-sarcomere at the start of the simulation, in nanometers.
*  pCa: Negative log of molar calcium content at the start of the simulation.
*  DataParams: Parameters relating to simulation characteristics rather than intended characteristics of the simulated muscle.
  *  DataSpan: The fraction of time that should be considered "steady state" for the purposes of data collection. For example, if .1 (the default) the last tenth of the simulated time will be considered steady state. Relevant variables should not fluctuate in this period.
  *  dt: Simulated time step, in seconds.
* knockout: Parameters relating to XB and RU knockout status.
  *  TnKOType: Boolean indicating the spatial type of knockout of RUs. 0 indicates "random" knockout, and 1 "uniform" along the filament free end.
  *  TnFraction: Fraction of active RUs. 1 means there is no knockout, 0 means there are no RUs at all, and .9 means a tenth of the RUs are knocked out.
  *  XBKOType: As TnKOType, but for cross-bridges.
  *  XBFraction: As TnFraction, but for cross-bridges.
* StiffScale: Scaling constants for the various springs.
  *  kxscaler: cross-bridge spring
  *  kascaler: actin (thin filament) spring
  *  kmscaler: myosin (thick filament) spring
*  filaments: pretty miscellaneous info.
  *  Tm\_Type: The type of RU cooperativity to use. See the DispatchCaRegCoop function.
  *  LACT, LMYO, L\_TITIN: Lengths of actin, myosin, and titin. For titin it's the length of the filament. I forget the other two (FIXME)
  *  NACT, NMYO: Number of thin and thick filaments (respectively)
  *  NumBridges: Bridges per filament. An extra 1 should be added as the end bridge.
  *  NumSites: Number of RU sites per filament. An extra 1 should be added as the end site.
  *  N\_Thick_Start: Number of myosin "heads" per XB ("crown").
  *  NTn: I forget (FIXME)
  *  KACT\_0, KMYO\_0: Spring constants of actin (thin) and myosin (thick), before scaling with StiffScale. In pN/nm.
  *  K\_TITIN: Spring constant of titin. I think in pN/nm? FIXME
  *  Angle\_Crowns: Pitch change between XBs, I think. In degrees.
  *  Angle\_Thick\_Start: Separation of heads on crowns, in degrees.
*  tcparam: Thermochemical parameters for XB state changes. These are used to derive the values of the 'thermochem' structure passed to OneStep.
  *  ATP, ADP, phos: Intracellular concentrations of ATP, ADP, and phosphate, in molar.
*  coop: Cooperativity information.
  *  Variables are named simply: XB2TF21 means the influence of XBs in state 2 on the RU transition from state 2 to state 1.
  *  In most cases these will be roots of a variable "Psi", because they are multiplicative. For example, by default XB2TF23, XB3TF23, and TF3TF23 are all the third root of Psi, so the TF 2-3 transition will be influenced by the third root of Psi cubed, which is of course 100.
  *  Coop\_Type: Form of cooperativity to use. This should always be 7, and in the future cooperativity changes can just be done by changing the variables already mentioned. For example, to have only RURU cooperativity, you should just set all the variables not in the form of "TF3..." to one.
*  TFRateScale: Scaling constants for RU state transitions. See ScaleThinFilRates and the CaRegCoop functions.

### Internal functions

#### OneRun

OneRun is the underlying workhorse of the system, containing the actual calls to produce one run. The parameters are essentially the same as those passed to the higher-level functions, but unpacked from structs, resulting in an ugly call convention with several dozen arguments.

#### DispatchCaRegCoop and OneStep

These are the slowest functions in the system, accounting for a majority of the runtime in most runs. DispatchCaRegCoop calls an underlying function based on the Tm_Type value (this underlying function is where most of the time is taken) that handles state changes in regulatory units, while OneStep handles state changes in the cross-bridges. These functions are run several thousand times for any run, once a simulated millisecond, and so even minor optimizations can have noticeable effects on runtime. For example, rewriting an if/else chain in a CaRegCoop function as a switch/case saved about ten seconds, even though these have the same semantics.

## Development

### Source control

Development has been moved to be under source control with the [git](http://git-scm.com/) system. This makes a history of changes as well as numerous other features available to developers. Git programs for Windows include [msys-git](https://msysgit.github.io/), [Github for Windows](https://windows.github.com/), and [SourceTree](http://www.sourcetreeapp.com/). Any one of these will work with the code or copies of the code as long as the .git and .gitignore folders and files are maintained. Editing of the code is still done through any text editor you like; git is a separate process that is run to "commit" already written changes.

Git has an instructive book written about it, [available free online](http://git-scm.com/book). The following very briefly explains a few basics, but for sustained use the book would really help.

In using source control, you can continue editing and running code through arbitrary editors (e.g., MATLAB's). Git is not an editor. Git is a separate program used when you want to manipulate and/or view the history of changes of the code.

The history can be visualized in msys-git by selecting "Repository -> Visualize All Branch History" in its menu. This view includes a graph of file changes over time, listings of the actual changes if one is selected, authorship information, etc.

You can "commit" a change, meaning it is added to the history. This is done in the Git interface by selecting changed files to "stage", and then "commit"ing them along with a brief message explaining the change. Files are selected and deselected using the functions in the "Commit" menu. "Rescan" tells msys-git to recheck the current folder for any changes you've made while msys-git was open.

Historical commits can be "checked out". This changes the files in the directory to be as they were at the time of that commit. This can be done temporarily, or permanently, in which case it is called a "reset". You can also "reset" specifically to the last committed version, obliterating any working changes. This feature obviously allows mistakes to be fixed easily without special backups.

History can be "branched". This allows multiple versions of the codebase to be worked on simultaneously. For example, you could have one branch adding one feature, and another removing an obsolescent one, and switch between them as you work so as to avoid having to think about both at once. Viewing an individual branch is also called "checking out" and works similarly to the historicla case. Branches are created and checked out with the "Branch" menu.

When branches are complete they can be "merged" back into one branch. This is done through the "Merge" menu. Watch out for merge conflicts which occur when two branches have changed the same file in incompatible ways.

### Future directions

#### Cleanup

Much of the source could stand to be reorganized. The CaRegCoop functions are full of redundant code that should be streamlined, and the Coop_Type parameter should be eliminated entirely as previously mentioned. The structs should be more documented and it should be decided whether unpacking them into variables is actually necessary.

#### Optimizations

As mentioned above, CaRegCoop functions and OneStep are the primary targets of optimization. Both should in principle be able to be run massively in parallel, as by GPGPU (general-purpose Graphics Processing Unit computing). This will require a reinterpretation of the semantics of these functions; e.g. see "occupied" in OneStep. MATLAB has GPGPU tools available, but they are limited to Nvidia's CUDA system. A better alternative would be the OpenCL cross-platform standard, but it is not available for MATLAB yet, so this probably would have to be done in C.

#### Data storage

Presently data from simulations is written into tabular text files with WriteText. Each parameter set is arbitrarily associated with a 'set number' making up a folder name to access. Improvements in this area are possible. Ideally, it would be possible to enter a parameter set and get previously produced data for those parameters, without having to worry about "set numbers". This could be done through some kind of persistent database format.

It might also be possible to use a more organized format than the present text files, but it's hard to beat the simplicity of tab-separated values.

#### Expanded possibilities

It should be possible to adjust the code to work with other models of state; that is, adding more RU or XB states and associated rate constants.