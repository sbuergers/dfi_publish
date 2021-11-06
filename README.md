# README

Experiment and analysis scripts accompanying the paper "The role of alpha oscillations in temporal binding within and across the senses" (Buergers & Noppeney, 2022).

DFI stands for "double flash illusion", although the experiment really investigated 1 vs 2 flash discrimination in various auditory context (0-2 beeps). 

## What is this repository for?

If you request the preprocessed data from this experiment you should be able to reproduce all analyses, results and figures in the manuscript. Unfortunately, we did not manage to find a storage solution where we could upload the data as used by the scripts. Instead, we low-pass filtered the data at 28 Hz and downsampled to 76 Hz. For the analyses, this should not make a huge difference, but it will break some of the scripts (at the very least due to different filenames) and maybe very slightly alter the results.

The code was not written with the intent of creating a toolbox or making it usable by a wider audience. I tried to clean it up enough so that it is understandable and relatively easy to follow which scripts belong to which analysis and the order in which they should be executed. 

Note that Figure and analysis scripts are sometimes entwined, if you do not find the analysis you are looking for in the seemingly designated folder (e.g. `eeg` for a particular eeg analysis), please look in the `fig` folder.

## Whom do I talk to?

Steffen Buergers
sbuergers@gmail.com
