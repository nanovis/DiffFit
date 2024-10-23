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

1. Open the source structure and the target map in ChimeraX. Let's use [`PDB-6WTI`](https://www.rcsb.org/structure/6wti) and its associated EM Map [`EMD-21897`](https://www.ebi.ac.uk/emdb/EMD-21897) as an example.
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
9. Save a molecule by clicking `Structure` if desired; and zero the density occupied by the current molecule:
    1. Click `Simulate volume`
    2. Change the surface level threshold for the simulated volume if necessary
    3. Click `Zero density`
7. Repeat the last step (Save, Simulate, Zero) for other plausible fits
8. Save the last `working volume` by `File > Save > Files of type: MRC > Map: working volume` 
   and use it as the input for another round if needed.
9. In the `Interactive` tab, you may click `Options` to change the fitting parameters.
10. In the `Settings` tab, you may change the global settings of DiffFit. The benchmark table in our paper uses `Fit atoms: All atoms` on an Nvidia RTX 4090 GPU.
11. You can find the computing time in the ChimeraX log window, usually at the right. Look for `DiffFit total time elapsed: `. The first run is slightly slower as there are some global initialization processes.
12. If you are interested in comparing DiffFit with the ChimeraX Fit in Map command, you may run a command similar to `fit #1 in #2 search 1000`. Check the [fitmap doc](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/fitmap.html) for more.



### Scenario 2: Composite multiple structures

The logic of compositing is you fit multiple structures in one run to one Cryo-EM volume map. 
Apart from the target volume file, you will also need to prepare all the individual structures
that you want to fit and simulate a map for each structure. 
You will run this functionality in `Disk` mode by specify the path to the involved files or folders.

As a demo, let's composite the individual chains from [PDB-8SMK](https://www.rcsb.org/structure/8SMK) 
into [EMD-40589](https://www.ebi.ac.uk/emdb/EMD-40589).   

1. Open `8SMK` in ChimeraX

2. Go to the `Utilities` tab. Under the `Split a structure into individual chains` section:
   1. Set an `Output Folder`. The default is a new folder called `split_out` in your working directory, 
      which usually is your desktop.
   2. Select `8smk.cif` as the `Structure`.
   3. Press `Split`.
   4. Optional: Go to the output folder, delete chain D, E, and F. Because they are the same as chain A, B, and C.

3. Go to the `Utilities` tab. Under the `Simulate a map for each structure in the folder` section:
   1. If your individual chains are in `split_out`. Then you may directly click `Simulate`. 
      By default, the simulated maps are in a new folder called `sim_out` in your working directory, 
      which usually is your desktop. 
   2. We explain all the parameters at the end. Please check them if you want to change the default values.

4. Go to the `Disk` tab, set all the file and folder paths, set the `Target Surface Threshold` 
   (you may open the map in the same ChimeraX window to decide a good surface level value).
5. Click `Run!`.
6. Use the same way as in `Scenario 1: Fit a single structure` to view the results.
7. Save the working volume and go for another round of fitting if needed. 
   

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
