## Reproduce guide

This guide is to reproduce Fig. 4 of our paper. 

---

### Video guide

### Text guide

1. Follow the [installation guide](https://github.com/nanovis/DiffFit?tab=readme-ov-file#install)
2. Download the [repository](https://github.com/nanovis/DiffFit/archive/refs/heads/main.zip) and unzip it to Desktop
3. Open ChimeraX, launch DiffFit, go to the `Disk` tab
4. Set the file path as below:
   1. Target Volume: `DiffFit-main\dev_data\input\domain_fit_demo_3domains\density2.mrc`
   2. Structures Folder: `DiffFit-main\dev_data\input\domain_fit_demo_3domains\subunits_cif`
   3. Structures Sim-map Folder: `DiffFit-main\dev_data\input\domain_fit_demo_3domains\subunits_mrc`
5. Click `Run!`
6. With a RTX 4090 GPU, it takes about 14 seconds to finish the computation. 
   If there is no CUDA-enabled GPU, the computation might take a few minutes to finish. 

---

Follow a similar process and with the help of the 
[Doc](https://github.com/nanovis/DiffFit/blob/main/Doc.md), 
one should also be able to roughly reproduce Fig. 1, 5, 6, 7, and Table 1, 2, 3. 

To get the results as close to our paper as possible: 
1. please use [DiffFit-0.6.0-py3-none-any.whl](https://github.com/nanovis/DiffFit/releases/download/v0.6.0/DiffFit-0.6.0-py3-none-any.whl)
2. go to the `Settings` tab, change `Fit atoms:` to `All atoms`. 