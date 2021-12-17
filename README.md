# Lymphomoid’s IF analysis pipeline
This pipeline was put together by Daniele Tavernari and Marco Varrone in December 2021 to streamline the analysis of multicolor immunofluorescence (IF) images of lymphoma and lymphomoid samples. It involves both manual and automated steps using open source tools and custom scripts that we developed. The pipeline is not optimal and can be improved in some of its parts. Feel free to adapt / modify / further develop the pipeline itself or the scripts involved, and reach out to us if you want to discuss.

## Step 1 - Image loading and inspection on QuPath
A good starting point is to load and visualize your IF image on QuPath. You can download QuPath and browse its documentation [here](https://qupath.github.io/). You can install QuPath on your local machine and open images stored on the server by mounting the corresponding folders on your local machine. Be aware that huge images require a good internet connection to be loaded and inspected without lags. A good amount of RAM is also advisable. 

QuPath is organized in “projects”. When you create a new project, you will be asked to allocate an empty folder to it. Then, you can add images to the project. It is fast to switch between images of the same project, and you will see the thumbnails of all of them at the same time; conversely, to inspect the images of a different project you will have to close the current one and open the other. This can affect the way you want to organize your project(s): I recommend to put related images into the same project (e.g. all images of the same patient, or of the same treatment). It is in principle also possible to store all images of all patients in the same project, but it probably would become messy.

When adding a new image, change its name into a meaningful one (sometimes the name of the raw file is just e.g. “Image_10.vsi” - you might want to change it into e.g. “HLS23_#23_multiplex”). 

## Step 2 - Lymphomoid(s) boundary drawing

## Step 3 - Fluorescence threshold optimization on the image 

## Step 4 - Nuclei detection and cell-level quantification with DeepCell

## Step 5 - (downstream analyses)

