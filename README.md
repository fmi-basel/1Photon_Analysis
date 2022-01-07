# 1Photon_Analysis

# Introduction #
This code is mainly a wrapper for several amazing projects and streamlines their usage to generate a coherent data output. If you want to learn more about the underlying packages, you can read about the [motion correction](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-33-2-156), [CNMFE](https://elifesciences.org/articles/28728) and the [Sheintuch](https://www.cell.com/cell-reports/pdf/S2211-1247(17)31430-4.pdf) and [Ahanonu](https://www.science.org/doi/full/10.1126/science.aap8586) cross-day registration. You can also check out the code, [here](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation), [here](https://github.com/zhoupc/CNMF_E), [here](https://github.com/zivlab/CellReg) and [here](https://github.com/bahanonu/ciatah).

For a complete, maintained implementation of many miniscope analysis tools, I highly recommend the package of Biafra Ahanonu - [CIAtah](https://github.com/bahanonu/ciatah). It comes with a wide range of GUIs that help with navigation for inexperienced users. 

This code repository serves as a minimal implementation (or parts) of the forenamed repositories and can be a good tool for people that look to customize their own pipelines or are looking for a code that is easy to run and works well out of the box without much parameter tuning. It has been tested on data from [MLR](https://www.cell.com/cell/pdf/S0092-8674(21)00828-X.pdf), [PFC](https://www.imb.de/students-postdocs/international-phd-programme/ipp-groups/beat-lutz), [AuC](https://www.physiologie.uni-freiburg.de/research-groups/ag-letzkus), [NAc](https://www.imb.de/students-postdocs/international-phd-programme/ipp-groups/beat-lutz) and [BLA](https://science.org/doi/10.1126/science.abg7277).  

The code currently works with the following data formats: tif, tiff, hdf5, mat and isxd (provided a current ISDP software installation exists).

# Installation #
The easiest way to use this code is to clone the repository in Matlab. You can check out how to do it [here](https://www.mathworks.com/help/matlab/matlab_prog/retrieve-from-git-repository.html) and [here](https://www.youtube.com/watch?v=O7A27uMduo0). Alternatively you can also download the code and save and maintain it locally on your computer. 

It requires a computer with a minimum of 64 GB RAM and for ease of processing should have at least 4 - 6 cores. It has been tested with Matlab 2017b and Matlab 2018a. It requires the following Matlab packages:  
  
* 'Optimization Toolbox'  
* 'Signal Processing Toolbox'  
* 'Image Processing Toolbox'  
* 'Statistics and Machine Learning Toolbox'  
* 'Curve Fitting Toolbox'  
* 'Parallel Computing Toolbox'  

To run the code, open the script *oneP_Image_Analysis.m* in your editor and follow the instructions below.

# Usage #
## Parameter Selection ##
Parameters are described extensively in inline comments in the code, but most importantly you need to identify the path of your Inscopix Installation when processing .isxd files from NV3 systems. All other parameters should work in the specified range and can be changed if results do not match expectation, but are generally suited to obtain good first results. Moreover, the code expects a numerical identifier in the folder names that is later used to save individual animals in different folders. Please specify the length of the identifier in the variable: length_identifier. If the code doesn't find an numerical identifier, it will save the last processing step in an undifferentiated folder, so when processing multiple animals, it will overwrite results.

## Session Selection and Motion Correction ##
After setting the parameters, running the code will open a prompt for every animal. This requires the user to add all sessions of a particular animal that should be processed. The code expects to receive folder directories that contain the imaging files of individual sessions (one folder/session). It treats all imaging files with the specified file extension in each folder as files originating from the same experiment and will sort them according to their recording timestamps. The order in which sessions are arranged will correspond to the order in the concatenation procedure later.  

Cross-day sessions should **always** be processed independently and not concatenated.

After session selection, another prompt asks the user for every animal to draw a rectangle onto one of the first images of the recording. This rectangle is the ROI that is used for motion correction and should be large, but contained within the visible boundries of the GRIN lens or optic window. Once you are satisfied with the selection, a double click will either open the prompt for the next animal, or let the code proceed to the motion correction stage.

First, the data is converted to a .mat file, which is used for later processing stages. After this step, the code will run through the video in small pieces and run the motion correction over the video as long as it takes to fall under the user specified threshold for motion in the video. Importantly the algorithm just deals with rigid motion and not with non-rigid motion, which can leave few videos with motion artifacts. From personal experience, optimization of surgical protocols and recording procedures is advisable over trying to rescue videos, but exceptions remain. In case you want to try rescuing videos from non-rigid motion I would recommend the package [NormCorre](https://github.com/flatironinstitute/NoRMCorre). Translational motion is calculated on filtered videos and applied to raw videos since the filtering introduces artifacts that can lead to easy misinterpretations later (most prominently sharp, negative baseline deflections).

After finishing the Motion Correction all videos will have a downsampled video visualisation and max intensity projections for visual inspection. Here it is important that neurons in the center appear round with clear edges and that, focusing on major landmark in the video (i.e. blood vessels), there is no visible translational movement left. You can click through all video's by clicking enter after clicking into the command window. If you are satisified with all video's you can click Yes on the next prompt and the code will proceed.
If you are not satisfied with the outcome, press No and restart the code from the beginning. Currently there is no way to process individual Sessions and you will need to re-run the entire batch of video's. The main issue for suboptimal Motion Correction is a ROI that is too small. Try increasing the ROI size and if it still doesn't work, check for visible non-rigid Motion or individual black frames which are most often the issue.

## CNMFE ##
The code will then run all raw, motion corrected video's through an unedited version of Pengcheng Zhou's [CNMFE](https://github.com/zhoupc/CNMF_E) implementation. Here an important feature is the selection of parallel pools that can speed up the process significantly, but also require large amounts of RAM, therefore a balance needs to be found. 

## CNMFE - Postselection ##
During post-Selection the user will have the chance to sort all components that are not automatically excluded manually with a GUI that is slightly edited from the CNMFE implementation. Follow the instruction displayed in the command window to proceed.

Important features to look out for are:

**Spatial component**  
(1) Components that have a roundish shape (paticularly in the center) and are not strongly elongated  
(2) Have clear borders and are not smeared out  
(3) Are not too small or too large in comparison with the other extracted components  
(4) Are not located within the region of the motion correction artifacts (visible as white bands)    
(5) Are located on top of a visible peak in the local correlation image

**Temporal component**  
(1) Have clear transients that conform with the biophysical properties of the calcium indicator (for example have a clearly visible decay that is consistent with the dissociation constant of the indicator)  
(2) Have a stable baseline that does not change abruptly (negative baseline changes are only biophysically feasible if they have a kinetic slower or equal to the dissociation constant of the calcium sensor, any faster transients must be artifact driven)  
(3) Have a good signal to noise ratio (SNR) - a line indicates mean + 3 * std, which can be a good first measure to assess SNR

After finishing all components of all sessions and animals the code will proceed to cross-day alignment.

## Cross-Day alignment ##
Here I used the code implementations of [Ahanonu](https://github.com/bahanonu/ciatah) and [Sheintuch](https://github.com/zivlab/CellReg), which will automatically align the components of all sessions that were selected and processed.

## Cross-Day Alignment verification ##
This step serves as a final sanity check to verify the automatic output of the Cross-Day alignment results and users are asked to manually accept all identified aligned components. In this step components that are only misaligned in individual sessions can also be excluded by indicating the Session in the command window. Features that we deem important for successfull alignment are the consistency in transient shape and spatial component overlap, both in spatial location as well as the appearance.

# Understanding Results #
The code will generate a new folder in each session that was processed containing all relevant intermediate processing steps and results. Folders contain the following variables:

`dataset.mat` - contains the raw, spatially downsampled data in form *x * y * t*  
`MC.mat` - contains the raw, spatially downsampled and motion corrected data in form *x * y * t*  
`MC_shifts.mat` - contains the applied pixel shifts during Motion Correction *(x, y) * t*  
`Template.mat` - Template that was used for motion correction (median projection image of the first 100 images after motion correction) *x * y*    
`CN.mat` - local pixel correlation Image, generated by the CNMFE *x * y*    
`Results_CNMFE.mat` - Results of the CNMFE extraction process (still in Sources2D format) *Source2D structure*   
Folder : `MC_Source_Extraction` - contains data relevant for CNMFE (no relevance for users)
`Results_CNMFE_postprocessed.mat` - Remaining components of CNMFE after manual selection and automatic exclusion (still in Sources2D format) *Source2D structure*    
`CNMFE_Data.mat` - subselected data from CNMFE after post-selection is finished, saved in normal structure (see below for list of variables) *structure*    
`Delete_1/2.mat` - Deletion indices for automatic exclusion *x*    
`Outline.mat/png` - Outline of the selected neurons (plotted on top of the correlation image)   

After the Cross-Day Alignment is finished a new folder is generated called *Alignment*. In this folder the code will save several variables related to the alignment procedure and the final Outputs *Data_Miniscope.mat* and *Data_Miniscope_PP.mat*. Both variables are identical and contain all information necessary for downstream analysis, however the variable **Data_Miniscope_PP.mat** is the Output generated, after the user finishes the post cross-day Alignment step.  It contains variables that follow the CNMFE naming convention set out by Pencheng Zhou. I recommend to follow his [advice](https://github.com/zhoupc/CNMF_E/wiki/Understand-CNMF-E-results) on utilizing the variables **C** and **A** for downstream analysis, however all variables are accessible.

To help with analysing cross-day data, the code will generate multiple variables containing the concatenated results in single variables:  

C - contains the variable C for components that could be aligned across all Sessions *neurons * t*  
C_raw - contains the variable C_raw for components that could be aligned across all Sessions *neurons * t*  
S - contains the variable S for components that could be aligned across all Sessions *neurons * t*  

C_all - contains the concatenated variable C for all global components. Sessions that could not be aligned are filled with nan's *neurons * t*  
C_raw_all - contains the concatenated variable C_Raw for all global components. Sessions that could not be aligned are filled with nan's *neurons * t*  
S_all - contains the concatenated variable S for all global components. Sessions that could not be aligned are filled with nan's *neurons * t*

C_unsorted - contains the variable C split into cells (one cell/session) with all neurons detected in this session *session{n}(neurons * t)*  
C_raw_unsorted - contains the variable C_raw split into cells (one cell/session) with all neurons detected in this session *session{n}(neurons * t)*  
S_unsorted - contains the variable S split into cells (one cell/session) with all neurons detected in this session *session{n}(neurons * t)*  

Order - contains the indices of neurons that could be aligned across days (example *Row*: 0 0 5 4 1 - this neuron could be aligned in 3 Sessions and can be called with C_unsorted{3}(5, :) for session 3, C_unsorted{4}(4, :) for session 4 and C_unsorted{5}(1, :) for Session 5). *neuron * session*  

Start_Frame_Session - contains the index of the start of each of the Session referring to the variables *..._all* and *C, C_raw and S*.  

SessionNames - contains the filepaths to the sessions that are concatenated (in order of concatenation).

# References #
If you decide to use this code, please cite the relevant works especially of the people who wrote the original Code that runs this pipeline:    

Zhou, P., Resendez, S.L., Rodriguez-Romaguera, J., Jimenez, J.C, Neufeld, S.Q., Giovannucci, A., Friedrich, J., Pnevmatikakis, E.A., Stuber, Garret D , Stuber, G.D., Hen, R., Kheirbek, M.A., Sabatini, B.L., Kass, R.E., Paninski, L. (2018). [Efficient and accurate extraction of in vivo calcium signals from microendoscopic video data.](https://elifesciences.org/articles/28728) eLife, pp.e28728.

Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, [Efficient subpixel image registration algorithms](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-33-2-156), Opt. Lett. 33, 156-158 (2008)

Sheintuch, L., Rubin, A., Brande-Eilat, N., Geva, N., Sadeh, N., Pinchasof, O., Ziv, Y. (2017). [Tracking the Same Neurons across Multiple Days in Ca2+ Imaging Data.](https://www.sciencedirect.com/science/article/pii/S2211124717314304) Cell Reports, 21(4), pp. 1102–1115. doi: 10.1016/j.celrep.2017.10.013.

Corder, G., Ahanonu, B., Grewe, B. F., Wang, D., Schnitzer, M. J., & Scherrer, G. (2019). [An amygdalar neural ensemble that encodes the unpleasantness of pain.](https://www.science.org/doi/full/10.1126/science.aap8586) Science, 363(6424), 276-281.

if you find this pipeline useful and use it in parts or whole, please cite also the following paper:

Courtin, J., Bitterman, Y., Müller, S., Hinz, J., Hagihara, K. M., Müller, C., Lüthi, A. (2022). [A neuronal mechanism for motivational control of behavior.](https://science.org/doi/10.1126/science.abg7277) Science, 375(6576)

# Help #
If you have further questions, please reach out! You can find my contact information [here](https://hijul.github.io/).

