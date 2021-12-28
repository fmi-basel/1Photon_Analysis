# 1Photon_Analysis

# Introduction #

This code is mainly a wrapper for several other amazing projects and streamlines it's usage to generate a coherent data output. If you want to learn more about the underlying packages you can read about the [Motion Correction](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-33-2-156), [CNMFE](https://elifesciences.org/articles/28728) and the [Sheintuch](https://www.cell.com/cell-reports/pdf/S2211-1247(17)31430-4.pdf) and [Ahanou](https://www.science.org/doi/full/10.1126/science.aap8586) cross-day registration. You can also check out the code, [here](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation), [here](https://github.com/zhoupc/CNMF_E), [here](https://github.com/zivlab/CellReg) and [here](https://github.com/bahanonu/ciatah).

For a complete, mantained implementation of many Miniscope Analysis tools I highly recommend the package of Biafra Ahanonu [CIAtah](https://github.com/bahanonu/ciatah). It comes with a wide range of GUIs that help with navigation for inexperienced users. 

This code repository serves as a minimal implementation of the forenamed repositories and can be a good tool for people that look to customize their own pipelines or are looking for a code that is streamlined and works, well out of the box without much parameter tuning and has been tested on data from vHPC, PFC, AuC, NAc, BLA and CeA. 

# Installation #
The easiest way to use this code is to clone the repository in Matlab. You can check out how to do it [here](https://www.mathworks.com/help/matlab/matlab_prog/retrieve-from-git-repository.html) and [here](https://www.youtube.com/watch?v=O7A27uMduo0). Alternatively you can also download the code and save and mantain it locally. 

It requires a computer with a minimum of 64 GB RAM and for ease of processing should have at least 4 - 6 cores. It requires following Matlab packages:

'Optimization Toolbox'
'Signal Processing Toolbox'
'Image Processing Toolbox'
'Statistics and Machine Learning Toolbox'
'Curve Fitting Toolbox'
'Parallel Computing Toolbox'
'MATLAB Parallel Server'

# Usage #

## Parameter Selection ##





