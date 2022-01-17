# Lymphomoid’s IF analysis pipeline
This pipeline was put together by Daniele Tavernari and Marco Varrone in December 2021 to streamline the analysis of multicolor immunofluorescence (IF) images of lymphoma and lymphomoid samples. It involves both manual and automated steps using open source tools and custom scripts that we developed. The pipeline is not optimal and can be improved in some of its parts. Feel free to adapt / modify / further develop the pipeline itself or the scripts involved, and reach out to us if you want to discuss.

## Step 1 - Image loading and inspection on QuPath
A good starting point is to load and visualize your IF image on QuPath. You can download QuPath and browse its documentation [here](https://qupath.github.io/). You can install QuPath on your local machine and open images stored on the server by mounting the corresponding folders on your local machine. Be aware that huge images require a good internet connection to be loaded and inspected without lags. A good amount of RAM is also advisable. 

QuPath is organized into “projects”. When you create a new project, you will be asked to allocate an empty folder to it. Then, you can add images to the project. It is fast to switch between images of the same project, and you will see the thumbnails of all of them at the same time; conversely, to inspect the images of a different project you will have to close the current one and open the other. This can affect the way you want to organize your project(s): I recommend to put related images into the same project (e.g. all images of the same patient, or of the same treatment). It is in principle also possible to store all images of all patients in the same project, but it probably would become messy.

When loading a .vsi raw image file, sometimes the same file will contain more than one image (different acquisitions of the same sample, or even different samples). Discard those that are unusable (e.g. out of focus, or with too little tissue and only fatty holes). Instead, for the good images, change their QuPath image name into a meaningful one (sometimes the name of the raw file is just e.g. “Image_10.vsi” - you might want to change it into e.g. “HLS23_s12_acq01”, which mentions the Human Lymphoma Sample number, the paraffin section, and the acquisition ID of that section).  It is advisable to rename them in a clear and consistent way.

Now you can explore the loaded image. If needed, set its type into "Fluorescence". You can change brightness, contrast and show/hide the different channels, corresponding to the different stained proteins. Sometimes it is useful to look at a single channel in grayscale mode. 

In the first batch of images, the channel-protein association was: 
* DAPI = Nuclei staining
* FITC = CD20 (B cells)
* CY5 = KI67 (proliferation)
* CFP = CD4 (CD4 T cells)
* RFP = CD8 (CD8 T cells)
* Alexa 594 = CD68 (macrophages)


## Step 2 - Lymphomoid(s) boundary drawing
In this step, you will draw and save the coordinates of the boundaries of all lymphomoids present in the image, so that they will be processed separately in the downstream analyses, and pieces of tissue scattered around the gel will be discarded.
1. For each lymphomoid, draw a closed boundary with QuPath polygon tool
2. Open _getBoundary.groovy_ script in QuPath script editor (Automate -> Show script editor)
3. In the script, set the absolute path to your desired output directory (_OutDir_) and the image name (_ImageName_)
4. For each boundary:
   * Select the boundary itself (double click - a selected boundary is displayed in yellow)
   * Set the lymphomoid name in the script (_LymphomoidName_). Be careful in this passage: if the same lymphomoid appears in different images, they should all have the same _LymphomoidName_. In this way, different acquisitions (encoded with different _ImageName_) will be treated as replicated measurements of the same lymphomoid in the downstream analyses
   * Run the script


## Step 3 - Fluorescence threshold calibration on the image
In this step, you will tune the thresholds for each channel to classify a cell as "positive" or not for each of the proteins. To do so, a script will automatically detect nuclei and classify cells according to the thresholds that you give it as input (for DAPI and Ki67, cells will be classified according to the fluorescence levels inside the nuclei; for all other markers, the fluorescence levels will be thresholded in the cytoplasmic regions). Tune the threshold of each channel until you are satisfied with the classification. Finally, the script stores all the tuned thresholds in a text file, that will be used in the downstream analyses. Here are the detailed (sub)steps:
1. Draw a rectangular region that seems to contain cells positive for each of the channels
2. Open _calibrateThresholds.groovy_ script in QuPath script editor
3. In the script, set the absolute path to your desired output directory (_OutDir_) and the image name (_ImageName_) in the _Input_ section. These have to be the same as for the boundary drawing (Step 2, point 3).
4. Tune the DAPI threshold first:
   * In the 'Brigthness & contrast' menu, turn off all channels except for DAPI, and set it to grayscale
   * Hover the cursor on the nuclei to have an idea on what might be the threshold
   * Assign a reasonable threshold guess to the _DAPI_thresh_ variable in the _Input_ section of the script
   * Run the script and inspect the results visually. You can view/hide the detections in QuPath with View->Show detections (Keyboard shortcut: D)
   * If nuclei have been under- (or over-) called, change the DAPI threshold and re-run the script until you are satisfied with the results
5. Tune the FITC threshold:
   * In the 'Brigthness & contrast' menu, show the FITC channel
   * Hover the cursor on the FITC stained cytoplasms to have an idea on what might be the threshold
   * Assign a reasonable threshold guess to the _FITC_thresh_ variable in the _Input_ section of the script
   * In the _Calibration_ section of the script, uncomment the line related to the FITC classification (_setCellIntensityClassifications("Cytoplasm: FITC mean", FITC_thresh_) by removing the groovy comment symbol ('//') 
   * Run the script and inspect the results visually. Cells positive for FITC will appear in red, negative cells will be blue.
   * Tune the FITC threshold and re-run the script until you are satisfied with the results
   * Once you are done, comment out again the line related to the FITC classification. Don't ever comment out the line related to the DAPI nuclei detection.
6. Repeat step 5 for all the other channels. Remember that at every tuning round only two lines in the _Calibration_ section should be uncommented: the DAPI detection and the intensity classification of your channel of interest
7. If you wish, repeat all the steps from 1 to 6 drawing a different rectangular region, to test whether your thresholds would stay the same or not
8. At each run, the script is saving and overwriting the output file with all the thresholds. Thus, your last run should be done with all the thresholds already tuned, as the final output file will be saved with those. If you want to run the script without saving the file, comment out the _Saving_ section

## Step 4 - Nuclei detection and cell-level quantification with DeepCell

## Step 5 - (downstream analyses)

