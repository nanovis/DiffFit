**DiffFit**: Visually-Guided **Diff**erentiable **Fit**ting of Molecule Structures to Cryo-EM Map

[![DiffFit_Teaser_30sec](https://github.com/nanovis/DiffFit/blob/51b03a45ada78340f949acb0953503fbc72cfb18/images/DiffFit_Teaser_30sec.gif?raw=true)](https://github.com/nanovis/DiffFit/blob/51b03a45ada78340f949acb0953503fbc72cfb18/images/DiffFit_Teaser_30sec.gif)

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

IEEE VIS 2024 Submission [arXiv preprint](https://arxiv.org/abs/2404.02465), [Video](https://youtu.be/dWcHDWT9_mw), [OSF repo](https://osf.io/5tx4q/)

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

Now, DiffFit should be fully installed. 
1. Launch it via `Tools > Volume Data > DiffFit`.
2. Right-click inside the DiffFit panel to access its `Help` page, and put it `In Favorites Menu`.
3. By default, the DiffFit panel floats above the ChimeraX window. You may right-click, check `Dockable Tool`, and move the panel around to dock it (suggest docking at the less preferred side, which usually is the left side). Then you may right-click and click `Save Tool Position`. 

## Demo usage scenarios

### Scenario 1: Fit a single structure

1. Open the source structure and the target map in ChimeraX. Let's use [`6WTI`](https://www.rcsb.org/structure/6wti) and its associated EM Map [`EMD-21897`](https://www.ebi.ac.uk/emdb/EMD-21897) as an example.
   1. Option 1, via ChimeraX command line. Run the following command in the command line at the very bottom of the ChimeraX windows. 
      1. `open 6WTI`
      2. `open 21897 from emdb`
   2. Option 2, via download the files and drag and drop.
      1. You can find the official RCSB webpage for a PDB ID by either Googling the ID or composing a URL: https://www.rcsb.org/structure/6WTI. Change `6WTI` to another ID if needed.
      2. `Download Files > PDBx/mmCIF Format`
      3. Click the first `EMDB`, after `EM Map EMD-ID` to access its EM map webpage
      4. `Download > 3D volume (map.gz)`; unzip the downloaded map file.
      5. Drag and drop both files (`6wti.cif` and `emd_21897.map`) into the ChimeraX window.
2. Change the iso-surface threshold level of the map.
      1. The Volume Viewer is usually located at the bottom-right of the ChimeraX window. Move the slide to change.
      2. DiffFit is very robust against this parameter, so very often, you don't have to change it. 
      3. To get the best performance, change the level to a value where you can see some secondary structures (alpha-helices or beta sheets). 
4. In the DiffFit panel, go to the `Interactive` tab. Click `Fit`.
5. After computing, DiffFit will automatically go to the `View` tab and select the top fit. You may change the threshold values. Usually, the defaults work fine. We explain all the parameters at the end.
6. You can now click on the rows to go through the fitting results. You may sort the table by a different metric by clicking the header (by default, the table is sorted by `Density`).
7. Once you find a plausible fit, you may use ChimeraX's Fit in Map tool (or a command similar to `fit #1 in #2`) to refine the placement.
8. If you have a ground truth structure to compare with, you may open that structure and use a command similar to `rmsd #1 to #3` to calculate the RMSD. Check the [RMSD doc](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/rmsd.html) for more. 
9. In the `Interactive` tab, you may click `Options` to change the fitting parameters.
10. In the `Settings` tab, you may change the global settings of DiffFit. The benchmark table in our paper uses `Fit atoms: All atoms` on an Nvidia RTX 4090 GPU.
11. You can find the computing time in the ChimeraX log window, usually at the right. Look for `DiffFit total time elapsed: `. The first run is slightly slower as there are some global initialization processes.
12. If you are interested in comparing DiffFit with the ChimeraX Fit in Map command, you may run a command similar to `fit #1 in #2 search 1000`. Check the [fitmap doc](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/fitmap.html) for more.



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
