## Reproduce guide

This guide is to reproduce Fig. 4 of our paper. 

---

### Video guide

1. [Install](https://youtu.be/JkAL7-T1U-U), and [enable CUDA](https://youtu.be/LCiVb9nA_Xo) 
2. [Demo Usage Scenario 3: Identify unknown densities](https://youtu.be/q4Au3DQ4vHU)

### Text guide

1. Download [DiffFit-0.6.2-py3-none-any.whl](https://github.com/nanovis/DiffFit/releases/download/v0.6.2/DiffFit-0.6.2-py3-none-any.whl)
2. Follow the [installation guide](https://github.com/nanovis/DiffFit?tab=readme-ov-file#install) 
3. Download the [repository](https://github.com/nanovis/DiffFit/archive/refs/heads/main.zip) and unzip it to Desktop
4. Open ChimeraX, launch DiffFit, go to the `Settings` tab, change `Fit atoms:` to `All atoms`. 
5. Go to the `Disk` tab, set the parameters as below:
   1. Target Volume: `DiffFit-main\dev_data\input\domain_fit_demo_3domains\density2.mrc`
   2. Structures Folder: `DiffFit-main\dev_data\input\domain_fit_demo_3domains\subunits_cif`
   3. Structures Sim-map Folder: `DiffFit-main\dev_data\input\domain_fit_demo_3domains\subunits_mrc`
   4. Target Surface Threshold: `1.44`
6. Click `Run!`
7. With a RTX 4090 GPU, it takes about 14 seconds to finish the computation. 
   If there is no CUDA-enabled GPU, the computation might take a few minutes to finish.
8. When the computation finishes, DiffFit's View tab will open. 
   The whole view should be very similar to Fig. 4 now. 
   To get an even closer image, do the following:
   1. In the DiffFit View tab's table section, click the first row and then the second row.
   2. At the right of ChimeraX, there is a panel named `Models`. 
      One row should give you `I7M317_D1.pdb`. Stay with this row.
   3. At the right of ChimeraX, there is a panel named `Volume Viewer`. 
      Under it, there is a gray square next to an eye, click the gray square,
      set `Alpha channel` to `125`.   
   4. Drag and dock the DiffFit window to the left of ChimeraX.
   5. Hold the left mouse button in the viewport in the middle and drag to rotate so that the view is the same as in the Fig. 4.   
   6. Change the background to white by clicking `White`, the 6-th last button in the top panel of ChimeraX.
   7. The numbers in the table will be slightly different. 
      Because there are some random generators in the computation.
      (By default, DiffFit randomly seeds 1000 initial placements for each structure 
      and optimizes them to better placements.) 

---

Following a similar process and with the help of the detailed 
[documentation](https://github.com/nanovis/DiffFit/blob/main/Doc.md), 
one should also be able to roughly reproduce Fig. 1, 5, 6, 7, and Table 1, 2, 3. 
