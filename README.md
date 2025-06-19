# Multispectral-NLOS-Imaging
Code that supplements the paper "Isolating Signals in Passive Non-Line-of-Sight Imaging using Spectral Content" accepted to the International Conference on Computational Photography (ICCP) 2023. This code can run all of the experiments in the main paper. All code is in MATLAB 2023a and the datasets are stored in it as well.

The code is separated into several folders. "experiment_runs" contains scripts that will run the experiments shown in the paper. "important_functions" contains code that performs the most essential portions of the multispectral imaging (unmixing and reconstructions). This is useful to parse through if interested in understanding the implementation of the methods.

## MATLAB toolboxs:
Curve Fitting Toolbox <br>
Image Processing Toolbox <br>
Statistics and Machine Learning Toolbox <br>

## Contact
Please contact the main author (Connor Hashemi) through email at connor.hashemi1995@gmail.com for any questions or comments.

## Citation
If using the code, please use the following citation: <br>
<br>
Connor Hashemi, Rafael Avelar, and James R. Leger. "Isolating Signals in Passive Non-Line-of-Sight Imaging using Spectral Content." *2023 IEEE International Conference on Computational Photography (ICCP)*. IEEE, 2023. \[Accepted\]

## Special Acknowledgements:
We use some code to translate the spectral filters to actual colored images. [Link](https://www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions) to the file exchange and cited as: 

Jeff Mather (2023). Spectral and XYZ Color Functions (https://www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions), MATLAB Central File Exchange. Retrieved June 16, 2023.
