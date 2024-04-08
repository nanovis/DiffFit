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
