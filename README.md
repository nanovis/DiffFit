**DiffFit**: Visually-Guided **Diff**erentiable **Fit**ting of Molecule Structures to Cryo-EM Map

IEEE VIS 2024 Submission [Video](https://youtu.be/dWcHDWT9_mw), [OSF repo](https://osf.io/5tx4q/)

## Install 

1. Download the repository and unzip to a path. The following guide will use `D:\GIT\DiffFitViewer` as this path. 
2. Run ChimeraX command `devel build D:\GIT\DiffFitViewer; devel install D:\GIT\DiffFitViewer`
3. Open the system command line shell, install PyTorch, Biopython, mrcfile, scikit-learn to ChimeraX's Python
    1. Find ChimeraX's Python, you may find [this guide](https://www.cgl.ucsf.edu/chimerax/docs/devel/ides_debugging_profiling.html) useful. The following commands will use `C:\Users\luod\AppData\Local\ChimeraX\bin\python.exe`.
    2. Install PyTorch via the following command or according to its [official doc](https://pytorch.org/get-started/locally/)
       ```
       C:\Users\luod\AppData\Local\ChimeraX\bin\python.exe -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
       ```
    3. Install Biopython, mrcfile, scikit-learn via the following command
       ```
       C:\Users\luod\AppData\Local\ChimeraX\bin\python.exe -m pip install biopython mrcfile scikit-learn
       ```

Now, DiffFit should be fully installed. Launch it via `Tools > Volume Data > DiffFit`

![image](https://github.com/nanovis/DiffFitViewer/assets/8460424/7c44d942-6c03-40bd-9791-5383214bafc1)

Right-click in the panel to access DiffFit's help page. 

![image](https://github.com/nanovis/DiffFitViewer/assets/8460424/d94283cd-4b4c-4d69-aa20-68e78f81e75f)

## Demo usage scenarios

### Scenario 1: Fit a single structure

1. Download [PDB-8JGF](https://www.rcsb.org/structure/8JGF) and [EMD-36232](https://www.ebi.ac.uk/emdb/EMD-36232) 
   1. note the resolution as `2.7`Å from the webpage
   2. extract the map
   3. put the files (`8jgf.cif` and `emd_36232.map`) under, for example, `D:\GIT\DiffFitViewer\run\input\8JGF` 
2. Drop both files into ChimeraX, 
   1. take a note for the `pixel` value from the log, which represents the grid spacing for this volume, which is `1.04` in this case
   2. move and rotate the molecule and then save it as `8JGF_transformed.cif`. This step is only for demo purpose and is not necessary for real use cases
3. Put `8JGF_transformed.cif` under `D:\GIT\DiffFitViewer\run\input\8JGF\subunits_cif`
4. Simulate a map for the molecule
   1. Create two folders, `subunits_mrc` and `subunits_npy`, under `D:\GIT\DiffFitViewer\run\input\8JGF\`
   2. Open a new ChimeraX session and run `runscript "D:\GIT\DiffFitViewer\src\convert2mrc_npy.py" "D:\GIT\DiffFitViewer\run\input\8JGF\subunits_cif" "D:\GIT\DiffFitViewer\run\input\8JGF\subunits_mrc" "D:\GIT\DiffFitViewer\run\input\8JGF\subunits_npy" 2.7 1.04`
5. Run DiffFit. Set the parameters as follows and hit `Run!`
   1. Target volume: `D:\GIT\DiffFitViewer\run\input\8JGF\emd_36232.map`
   2. Structures folder: `D:\GIT\DiffFitViewer\run\input\8JGF\subunits_cif`
   3. Structures sim-map folder: `D:\GIT\DiffFitViewer\run\input\8JGF\subunits_mrc`
   4. Output folder: `D:\GIT\DiffFitViewer\run\output\8JGF`
   5. Experiment name: `fit_single_demo`
   6. Target surface threshold: `0.20`. Or use the author recommended contour level `0.162`. DiffFit is very robust against this parameter, a value between 0.02 - 0.4 is fine in this case.
   7. Leave the rest as default and hit `Run!`
6. After freezing for a few seconds (less than 15 seconds on RTX 4090), ChimeraX should be back and responsive to you. Click the `View` tab to examine the results.
   1. Save the molecule if desired
   2. You may take a look at the optimization steps
7. If you want to change the cluster tolerance, or if you run Compute on a cluster, or if you accidentally close ChimeraX after _Compute_ run, you can _View_ the results by the following parameter settings
   1. Target volume: `D:\GIT\DiffFitViewer\run\input\8JGF\emd_36232.map`
   2. Structures folder: `D:\GIT\DiffFitViewer\run\input\8JGF\subunits_cif`
   3. Data folder: `D:\GIT\DiffFitViewer\run\output\8JGF\fit_single_demo`
   4. Clustering - Shift Tolerance: `3.0` or the value you desire
   5. Clustering - Angle Tolerance: `6.0` or the value you desire
   6. Hit `Load`


### Scenario 2: Composite multiple structures


### Scenario 3: Identify unknown densities

