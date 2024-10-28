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

## YouTube tutorial videos

1. Install
   1. For users (demo version: v0.6.0): https://youtu.be/JkAL7-T1U-U
   2. For developers (demo version: before v0.2.0): https://youtu.be/aYqNZ0SNUfk
2. Demo Usage [Scenario 1: Fit a single structure](https://youtu.be/dHquT2Lsh54)
3. Demo Usage Scenario 2: Composite multiple structures
4. Demo Usage [Scenario 3: Identify unknown densities](https://youtu.be/4fV-qHO9spw)


## Install 

### Option 1 - From the official ChimeraX Toolshed: 
1. Download, install, and open [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html)
2. From the top menu bar, click `Tools > More Tools...`
3. From the newly opened ChimeraX's built-in browser page, find and click "DiffFit". You may click "more newest releases »" if you cannot find it on the home page.
4. Click `Install` (a big blue button). If it's your first time installing, it might take several minutes because it needs to install dependent packages such as PyTorch. 

### Option 2 - From the GitHub release page: 
1. Download the [latest distribution](https://github.com/nanovis/DiffFit/releases/latest) (download the `.whl` file)
2. Open ChimeraX and run the command `toolshed install <path to the downloaded .whl file>`

Now, DiffFit should be fully installed. 
1. Launch it via `Tools > Volume Data > DiffFit`.
2. Right-click inside the DiffFit panel to access its `Help` page, and put it `In Favorites Menu`.
3. By default, the DiffFit panel floats above the ChimeraX window. You may right-click, check `Dockable Tool`, and move the panel around to dock it (suggest docking at the less preferred side, which usually is the left side). Then you may right-click and click `Save Tool Position`. 

## Doc

Open the structure and the volume map, click `Fit`. 

---

The whole UI design is with the hope to allow users to be able to use DiffFit without training.
This is especially the case if the user is familiar with ChimeraX's 
[Fit in Map tool](https://www.cgl.ucsf.edu/chimerax/docs/user/tools/fitmap.html) 
or [fitmap command](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/fitmap.html).

We do also provide a detailed 
[Doc](https://github.com/nanovis/DiffFit/blob/main/Doc.md).

If you want to reproduce the figures in our paper, please refer to 
[Reproduce](https://github.com/nanovis/DiffFit/blob/main/Reproduce.md). 

If you see anywhere we can improve the design or the doc, please raise an 
[issue](https://github.com/nanovis/DiffFit/issues).  

