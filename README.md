**DiffFit**: Visually-Guided **Diff**erentiable **Fit**ting of Molecule Structures to Cryo-EM Map

IEEE VIS 2024 Submission [arXiv preprint](https://arxiv.org/abs/2404.02465), [Video](https://youtu.be/dWcHDWT9_mw), [OSF repo](https://osf.io/5tx4q/)

[![Image from Gyazo](https://i.gyazo.com/54a2b17ff4f914e6c87e78c01d4982dc.gif)](https://gyazo.com/54a2b17ff4f914e6c87e78c01d4982dc)

[![Image from Gyazo](https://i.gyazo.com/3bacb110ddbadb4cfffa34fb061b8f34.gif)](https://gyazo.com/3bacb110ddbadb4cfffa34fb061b8f34)

If you use material from this repository, please cite the associated paper:

Deng Luo, Zainab Alsuwaykit, Dawar Khan, Ondřej Strnad, Tobias Isenberg, and Ivan Viola. DiffFit: Visually-Guided Differentiable Fitting of Molecule Structures to a Cryo-EM Map. IEEE Transactions on Visualization and Computer Graphics, 31, 2025. To appear. doi: [10.1109/TVCG.2024.3456404](https://doi.org/10.1109/TVCG.2024.3456404)

bibTeX:
```
@article{Luo:2025:DVG,
  author      = {Deng Luo and Zainab Alsuwaykit and Dawar Khan and Ond{\v{r}}ej Strnad and Tobias Isenberg and Ivan Viola},
  title       = {{DiffFit}: Visually-Guided Differentiable Fitting of Molecule Structures to a Cryo-{EM} Map},
  journal     = {IEEE Transactions on Visualization and Computer Graphics},
  year        = {2025},
  volume      = {31},
  doi         = {10.1109/TVCG.2024.3456404},
  doi_url     = {https://doi.org/10.1109/TVCG.2024.3456404},
  github_url  = {https://github.com/nanovis/DiffFit},
  osf_url     = {https://osf.io/5tx4q/},
  preprint    = {https://doi.org/10.48550/arXiv.2404.02465},
  hal_url     = {https://hal.science/hal-04665408},
  video       = {https://youtu.be/dWcHDWT9_mw},
}
```

## YouTube tutorial videos (coming soon)

1. [Install](https://youtu.be/aYqNZ0SNUfk)
2. Demo Usage [Scenario 1: Fit a single structure](https://youtu.be/dHquT2Lsh54)
3. Demo Usage Scenario 2: Composite multiple structures
4. Demo Usage [Scenario 3: Identify unknown densities](https://youtu.be/4fV-qHO9spw)


## Install 

Option 1 - From the official ChimeraX Toolshed: 
1. Download, install, and open [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html)
2. From the top menu bar, click `Tools > More Tools...`
3. From the newly opened ChimeraX's built-in browser page, find and click "DiffFit". You may click "more newest releases »" if you cannot find it on the home page.
4. Click `Install` (a big blue button). If it's your first time installing, it might take several minutes because it needs to install dependent packages such as PyTorch. 

Option 2 - From the GitHub release page: 
1. Download the [latest distribution](https://github.com/nanovis/DiffFit/releases/latest) (download the `.whl` file)
2. Open ChimeraX and run the command `toolshed install <path to the downloaded .whl file>`

Now, DiffFit should be fully installed. Launch it via `Tools > Volume Data > DiffFit`. 

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
   2. move and rotate the molecule and then save it (select it, choose "Save selected atoms only", uncheck "Use untransformed coordinates") as `8JGF_transformed.cif`. This step is only for demo purpose and is not necessary for real use cases
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
6. After freezing for a couple of seconds (less than 15 seconds on one RTX 4090), ChimeraX should be back and responsive to you. Click the `View` tab to examine the results.
   1. Save the molecule if desired
   2. You may take a look at the optimization steps
7. If you want to change the cluster tolerance, or if you run Compute on a cluster, or if you accidentally close ChimeraX after _Compute_ run, you can _View_ the results by the following parameter settings
   1. Target volume: `D:\GIT\DiffFitViewer\run\input\8JGF\emd_36232.map`
   2. Structures folder: `D:\GIT\DiffFitViewer\run\input\8JGF\subunits_cif`
   3. Data folder: `D:\GIT\DiffFitViewer\run\output\8JGF\fit_single_demo`
   4. Clustering - Shift Tolerance: `0.5` or the value you desire
   5. Clustering - Angle Tolerance: `0.5` or the value you desire
   6. Hit `Load`


### Scenario 2: Composite multiple structures

1. Download [PDB-8SMK](https://www.rcsb.org/structure/8SMK) and [EMD-40589](https://www.ebi.ac.uk/emdb/EMD-40589) 
   1. note the resolution as `3.5`Å from the webpage
   2. extract the map
   3. put the files (`8smk.cif` and `emd_40589.map`) under, for example, `D:\GIT\DiffFitViewer\run\input\8SMK` 
2. Drop both files into ChimeraX, 
   1. take a note for the `pixel` value from the log, which represents the grid spacing for this volume, which is `0.835` in this case
   2. move and rotate the molecule and then save it (select it, choose "Save selected atoms only", uncheck "Use untransformed coordinates") as `8SMK_transformed.cif`. This step is only for demo purpose and is not necessary for real use cases
3. Create a folder `subunits` under `D:\GIT\DiffFitViewer\run\input\8SMK`
4. Split the chains into individual .cif files and simulate a map for each chain
   1. Open a new ChimeraX session and run `runscript "D:\GIT\DiffFitViewer\src\split_chains.py" "D:\GIT\DiffFitViewer\run\input\8SMK\8SMK_transformed.cif" "D:\GIT\DiffFitViewer\run\input\8SMK\subunits" 3.5 0.835`
   2. Put all generated .cif files under `D:\GIT\DiffFitViewer\run\input\8SMK\subunits_cif`
   3. Put all generated .mrc files under `D:\GIT\DiffFitViewer\run\input\8SMK\subunits_mrc`
   4. Delete all generated .npy files, or put them under `D:\GIT\DiffFitViewer\run\input\8SMK\subunits_npy`
   5. Keep only the unique chains (A, B, C) in `subunits_cif` and `subunits_mrc`
5. Run DiffFit. Set the parameters as follows and hit `Run!`
   1. Target volume: `D:\GIT\DiffFitViewer\run\input\8SMK\emd_40589.map`
   2. Structures folder: `D:\GIT\DiffFitViewer\run\input\8SMK\subunits_cif`
   3. Structures sim-map folder: `D:\GIT\DiffFitViewer\run\input\8SMK\subunits_mrc`
   4. Output folder: `D:\GIT\DiffFitViewer\run\output\8SMK`
   5. Experiment name: `round1`
   6. Target surface threshold: `0.8`. Or use the author recommended contour level `5.0`. DiffFit is very robust against this parameter, a value between 0.1 - 5.0 is fine in this case.
   7. \# shifts: `30`
   8. \# quaternions: `300`
   9. Leave the rest as default and hit `Run!`
6. After freezing for a couple of seconds (less than 30 seconds on one RTX 4090), ChimeraX should be back and responsive to you. Click the `View` tab to examine the results.
   1. Examine the fit, sort by a different metric
   2. If you want to change the cluster tolerance, or if you run Compute on a cluster, or if you accidentally close ChimeraX after _Compute_ run, you can _View_ the results by the following parameter settings
      1. Target volume: `D:\GIT\DiffFitViewer\run\input\8SMK\emd_40589.map`
      2. Structures folder: `D:\GIT\DiffFitViewer\run\input\8SMK\subunits_cif`
      3. Data folder: `D:\GIT\DiffFitViewer\run\output\8SMK\composite_unique_chains`
      4. Clustering - Shift Tolerance: `6` or the value you desire
      5. Clustering - Angle Tolerance: `15` or the value you desire
      6. Hit `Load`
   3. Save a molecule if desired
   4. Set the Resolution as `3.5`, and click `Simulate volume`
   5. Change the surface level threshold for the simulated volume if necessary
   6. Click `Zero density`
   7. Repeat the last 4 steps (Save, Simulate, Zero) for the same `Mol Id` at a different place, or for a different `Mol Id` until there is no good fit
   8. Save the last `working volume` by `File > Save > Files of type as MRC > Map as the desired one` as a new name, for example, `emd_40589_round_1.mrc`  
7. Repeat Step 5-6 until satisfied with the whole compositing
   1. Change the Target volume as: `D:\GIT\DiffFitViewer\run\input\8SMK\emd_40589_round_1.mrc`
   2. If needed, take out the already fitted chains from `subunits_cif` and `subunits_mrc`
   3. Give a new Experiment name: `round2`
   4. You may lower the \# shifts, for example, to `10`, and the \# quaternions to `100`
   5. Hit `Run!`
   

### Scenario 3: Identify unknown densities

The whole procedure is the same as in [Scenario 1: Fit a single structure](https://github.com/nanovis/DiffFitViewer?tab=readme-ov-file#scenario-1-fit-a-single-structure), 
only that there will be multiple structures under `subunits_cif`. 

There is a demo data set with one volume map and three structures to search against.
If you have put DiffFit under `D:\GIT\DiffFitViewer`, 
you can just hit `Run!` in the _Compute_ tab and then go to the _View_ tab. 
If otherwise, you just need to change the path for the input and the output data.

If you want to search against the whole candidate library for this case from DomainFit, 
you can either follow Steps 1-3 from its [doc](https://github.com/builab/DomainFit/blob/main/example/example.md)
to generate the PDB files for the domains, or just download the ones generated by us from 
[this Google Drive link](https://drive.google.com/file/d/1YT8puA1KBjnT9gWyt64-FWSy64wVaou0/view?usp=sharing). Of note is that we generated 359 PDB files by following DomainFit's Steps 1-3, 
instead of the mentioned 344 files. 

The computing time for searching the whole candidate library on one RTX 4090 is about 10 minutes. 
