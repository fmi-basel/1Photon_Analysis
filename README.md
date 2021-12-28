# 1Photon_Analysis

# Introduction #
This code is mainly a wrapper for several other amazing projects and streamlines it's usage to generate a coherent data output. If you want to learn more about the underlying packages you can read about the [Motion Correction](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-33-2-156), [CNMFE](https://elifesciences.org/articles/28728) and the [Sheintuch](https://www.cell.com/cell-reports/pdf/S2211-1247(17)31430-4.pdf) and [Ahanonu](https://www.science.org/doi/full/10.1126/science.aap8586) cross-day registration. You can also check out the code, [here](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation), [here](https://github.com/zhoupc/CNMF_E), [here](https://github.com/zivlab/CellReg) and [here](https://github.com/bahanonu/ciatah).

For a complete, mantained implementation of many Miniscope Analysis tools I highly recommend the package of Biafra Ahanonu [CIAtah](https://github.com/bahanonu/ciatah). It comes with a wide range of GUIs that help with navigation for inexperienced users. 

This code repository serves as a minimal implementation of the forenamed repositories and can be a good tool for people that look to customize their own pipelines or are looking for a code that is easy to run headless and works well out of the box without much parameter tuning. It has been tested on data from MLR, vHPC, PFC, AuC, NAc, BLA and CeA. 

# Installation #
The easiest way to use this code is to clone the repository in Matlab. You can check out how to do it [here](https://www.mathworks.com/help/matlab/matlab_prog/retrieve-from-git-repository.html) and [here](https://www.youtube.com/watch?v=O7A27uMduo0). Alternatively you can also download the code and save and mantain it locally on your computer. 

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
Parameters are described extensively in inline comments in the code, but most importantly you need to identify the path of your Inscopix Installation when processing .isxd files from NV3 systems. All other parameters should work in the specified range and can be changed if results do not match expectation, but are generally suited to obtain good first results.

## Session Selection and Motion Correction ##
After running the code, a prompt will appear for every animal that requires the user to add all sessions of a paticular animal that should be processed. The code expects to receive folder directories that contain the imaging files of individual sessions (One Folder/Session). It treats all imaging files with the specified ending in each folder as files originating from the same experiment and will sort them according to their recording timestamps. Cross-day sessions should **always** be processed independently and not concatenated.

After session selection another prompt asks the user for every animal to draw a rectangle onto one of the first Images of the recording. This rectangle is the ROI that is used for Motion Correction and should be large, but contained within the visible boundries of the GRIN lens or optic window. Once you are satisfied with the selection, a double click will either open the prompt for the second animal or let the code proceed to the Motion Correction Stage.

First, the data is converted to a .mat file, which is used for later processing stages. After this step the code will run through the video in chunks and run the motion correction over the video as long as it takes to fall under the user specified threshold for motion in the video. Importantly the algorithm just deals with rigid motion and not with non-rigid motion, which can leave few video's with severe Motion artifacts. From personal experience optimization of surgical protocols and recording procedures is advisable over trying to rescue video's, but exceptions remain. In case you want to try rescuing video's from non-rigid motion I would recommend the package [NormCorre](https://github.com/flatironinstitute/NoRMCorre). Also, shifts are calculated on filtered video's and applied to raw video's since the filtering introduces artifacts that can lead to easy misinterpretations later (most prominently sharp, negative baseline deflections).

After finishing the Motion Correction all videos will have a downsampled video visualisation and max intensity projections for visual inspection. Here it is important that neurons in the center appear round with clear edges and that, focusing on major landmark in the video (i.e. blood vessels), there is no visible translational movement left. You can click through all video's by clicking enter after clicking into the command window. If you are satisified with all video's you can click Yes on the next prompt and the code will proceed.
If you are not satisfied with the outcome, press No and restart the code from the beginning. Currently there is no way to process individual Sessions and you will need to re-run the entire batch of video's. The main issue for suboptimal Motion Correction is a ROI that is too small. Try increasing the ROI size and if it still doesn't work, check for visible non-rigig Motion or individual black frames which are most often the issue.

## CNMFE ##
The code will then run all raw, motion corrected video's through an unedited version of Pengcheng Zhou's [CNMFE](https://github.com/zhoupc/CNMF_E) implementation. Here an important feature is the selection of parallel pools that can speed up the process significantly, but also require large amounts of RAM, therefore a balance needs to be found. 

## CNMFE - Postselection ##
During post-Selection the user will have te chance to sort all components that are not automatically excluded manually with a GUI that is slightly edited from the CNMFE implementation. Follow the instruction displayed in the command window to proceed.

Important features to look out for are:

**Spatial Component**
(1) Components that have a roundish shape (paticularly in the center) and are not strongly elongated 
(2) Have clear borders and are not smeared out
(3) Are not too small or too large in comparison with the other extracted components
(4) Are not located within the region of the Motion Correction artifacts (visible as white bands)

**Temporal component**
(1) Have clear transients that conform with the biophsical properties of the calcium indicator (for example have a clearly visible decay that is consistent with the dissociation constant of the Indicator)
(2) Have a stable baseline that does not change abruptly (negative baseline changes can only biologically exist with the same dissociation constant as calcium transients)
(3) Have a good S/N ration (a line indicates mean + 3 * std, which can be a good first measure to assess S/N ratio)

After finishing all components of all sessions and animals the code will proceed to cross-day alignment.

## Cross-Day alignment ##
Here I used the code implementations of [Ahanonu](https://github.com/bahanonu/ciatah) and [Sheintuch](https://www.cell.com/cell-reports/pdf/S2211-1247(17)31430-4.pdf), which will automatically align the components of all Sessions that were selected and processed.

## Cross-Day Alignment verification ##
This step serves as a final sanity check to verify the automatic output of the Cross-Day alignment results and users are asked to manually accept all identified aligned components. In this step misaligned components can also be excluded by indicating the Session in the command window. features that we deem important for successfull alignment are the consistency in transient shape and spatial component overlap, both in spatial location as well as the appearance.

# Understanding Results #








