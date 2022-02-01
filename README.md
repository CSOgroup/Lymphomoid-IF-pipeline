# Lymphomoid’s IF analysis pipeline
This pipeline was put together by Daniele Tavernari and Marco Varrone in December 2021 to streamline the analysis of multicolor immunofluorescence (IF) images of lymphoma and lymphomoid samples. It involves both manual and automated steps using open source tools and custom scripts that we developed. The pipeline is not optimal and can be improved in some of its parts. Feel free to adapt / modify / further develop the pipeline itself or the scripts involved, and reach out to us if you want to discuss.

## Step 1 - Image loading and inspection on QuPath
A good starting point is to load and visualize your IF image on QuPath. You can download QuPath and browse its documentation [here](https://qupath.github.io/). You can install QuPath on your local machine and open images stored on the server by mounting the corresponding folders on your local machine. Be aware that huge images require a good internet connection to be loaded and inspected without lags. A good amount of RAM is also advisable. 

QuPath is organized into “projects”. When you create a new project, you will be asked to allocate an empty folder to it. Then, you can add images to the project. It is fast to switch between images of the same project, and you will see the thumbnails of all of them at the same time; conversely, to inspect the images of a different project you will have to close the current one and open the other. This can affect the way you want to organize your project(s): I recommend to put related images into the same project (e.g. all images of the same patient, or of the same treatment). It is in principle also possible to store all images of all patients in the same project, but it probably would become messy.

When loading a .vsi raw image file, sometimes the same file will contain more than one image (different acquisitions of the same sample, or even different samples). Discard those that are unusable (e.g. out of focus, or with too little tissue and only fatty holes). Instead, for the good images, change their QuPath image name into a meaningful one (sometimes the name of the raw file is just e.g. “Image_10.vsi” - you might want to change it into e.g. “HLS23_s12_acq01”, which mentions the Human Lymphoma Sample number, the paraffin section, and the acquisition ID of that section).  It is advisable to rename them in a clear and consistent way.

Now you can explore the loaded image. If needed, set its type into "Fluorescence". You can change brightness, contrast and show/hide the different channels, corresponding to the different stained proteins. Sometimes it is useful to look at a single channel in grayscale mode. 


## Step 2 - Lymphomoid(s) boundary drawing
In this step, you will draw and save the coordinates of the boundaries of all lymphomoids present in the image, so that they will be processed separately in the downstream analyses, and pieces of tissue scattered around the gel will be discarded.
1. For each lymphomoid, draw a closed boundary with QuPath polygon tool
2. Open _getBoundary.groovy_ script in QuPath script editor (Automate -> Show script editor)
3. In the script, set the absolute path to your desired main directory for all pipeline input/output files (_MainDir_) and the image name (_ImageName_)
4. For each boundary:
   * Select the boundary itself (double click - a selected boundary is displayed in yellow)
   * Set the patient-specific lymphomoid name in the script (_PatientLymphomoidName_), avoid using underscores. This name should include the patient ID and the lymphomoid ID (e.g. "LP02ctrl01"). Be careful in this passage: if the same lymphomoid appears in different images, they should all have the same _PatientLymphomoidName_. In this way, different acquisitions (encoded with different _ImageName_'s) will be treated as replicated measurements of the same lymphomoid in the downstream analyses
   * Run the script


## Step 3 - Fluorescence threshold calibration on the image
In this step, you will tune the thresholds for each channel to classify a cell as "positive" or not for each of the proteins. To do so, a script will automatically detect nuclei and classify cells according to the thresholds that you give it as input (for DAPI and Ki67, cells will be classified according to the fluorescence levels inside the nuclei; for all other markers, the fluorescence levels will be thresholded in the cytoplasmic regions). The script also takes as input a configuration table that stores the association between channels, antibodies and cellular location (e.g. FITC - B220 - Cytoplasm); the two provided configuration tables _human_channels.txt_ and _mouse_channels.txt_ should work in the majority of the cases, but check that they are correct before running the script and create new ones with the same structure in case it is needed. Then, tune the threshold of each channel until you are satisfied with the classification. Finally, the script stores all the tuned thresholds in a text file, that will be used in the downstream analyses. Here are the detailed (sub)steps:
1. Draw a rectangular region that seems to contain cells positive for each of the channels, and select it (you can also use an existing lymphomoid boundary)
2. Open _calibrateThresholds.groovy_ script in QuPath script editor
3. In the script, set the absolute path to your desired main directory (_MainDir_) and configuration table (_ConfigTable_), and the image name (_ImageName_) in the _Input_ section. These have to be the same as for the boundary drawing (Step 2, point 3).
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
   * In the _Calibration_ section of the script, uncomment the lines related to the FITC classification (_setCellIntensityClassifications(channel_celllocation_map['FITC']+(...)_) by removing the groovy comment symbol ('//') 
   * Run the script and inspect the results visually. Cells positive for FITC will appear in red, negative cells will be blue.
   * Tune the FITC threshold and re-run the script until you are satisfied with the results
   * Once you are done, comment out again the lines related to the FITC classification. Don't ever comment out the line related to the DAPI nuclei detection.
6. Repeat step 5 for all the other channels. Remember that at every tuning round only the lines related to (1) the DAPI detection and (2) the intensity classification of your channel of interest should be uncommented in the _Calibration_ section
7. If you wish, repeat all the steps from 1 to 6 drawing a different rectangular region, to test whether your thresholds would stay the same or not
8. At each run, the script is saving and overwriting the output file with all the thresholds. Thus, your last run should be done with all the thresholds already tuned, as the final output file will be saved with those. If you want to run the script without saving the file, comment out the _Saving_ section

## Step 4 - Nuclei detection and cell-level quantification with DeepCell (work in progress)
DeepCell is a nuclei and cell segmentation software that is more robust to different levels of marker intensity and, thus, gives better results when the intensity of DAPI varies dramatically in the same sample.

The whole process of detecting cells and extract the intensity of the markers is fairly automatized. The crucial step is to setup the input directories and files properly.

The script in run in Step 2 should have already created, inside the main directory, a directory called _Quantification_. _Quantification_ contains other directories with the names of the images for which boundaries have been drawn and finally, each of this directory should contain an empty directory called _registration_.

To summarize, you should find a structure like the following:

```bash
Quantification
├── HLS07_s01_acq01
│   ├── registration
├── HLS25_s41_acq01
│   ├── registration
├── HLS25_s41_acq03
│   ├── registration
...
```

### VSI to .OME.TIF format conversion
In this step we will convert the images from the .vsi to the .ome.tif format. The script uses a software called _bfconvert_ from Bio-Formats Tools.

1. If it doesn't already exist, create a directory called _Tiff_ inside the main directory. So it should be at the same level of the _Quantification_ directory.
2. Connect to the server using ssh on your Terminal and move to the directory containing the scripts of the pipeline.
3. Run the script _vsi2tiff.py_ with the following required parameters:
   * `--input_vsi`: the path to the .vsi file that needs to be converted into .ome.tif files.
   * `--output_dir`: the path to the _Tiff_ directory

   <details>
   <summary>Some other parameters may be required in case of non-standard uses of the pipeline, so they can usually be ignored</summary>

   * `--bftools_dir`: the path to the _bftools_ director. The default is `/mnt/data2/shared/Lymphomoid-IF-software/bftools` where it's already present.

   </details>
4. Note that a single vsi will generate a separate .ome.tif file for each image acquisition. The script will create multiple .ome.tif files concatenating the name of the parent directory of the .vsi file, the name of the .vsi file and the acquisition number.
5. Delete the acquisitions that have been discarded in Step 1. You can identify those acquisitions because there is no equivalent directory inside _Quantification_.
6. Rename and move each image into the _registration_ folder of the corresponding directory. So, every registration folder should contain one and only one .ome.tif file, as in the following diagram:
   ```bash
   Quantification
   ├── HLS07_s01_acq01
   │   ├── registration
   │   │   HLS07_s01_acq01.ome.tif
   ├── HLS25_s41_acq01
   │   ├── registration
   │   │   HLS25_s41_acq01.ome.tif
   ├── HLS25_s41_acq03
   │   ├── registration
   │   │   HLS25_s41_acq03.ome.tif
   ...
   ```

### Cell detection (i.e. segmentation)
The cell detection script require few Python packages for basic file processing. They are all already installed as a Python virtual environment in _/mnt/data2/shared/Lymphomoid-IF-software/Lymphomoid-IF-venv/_.
1. To be able to access the packages you need to activate the environment using:

   `source /mnt/data2/shared/Lymphomoid-IF-software/Lymphomoid-IF-venv/bin/activate`

2. Then you need to run the _cellDetection.py_ script, with the following required parameters:
   * `--sample_names`: the name of the images for cell detection, it could be 1 or many, separated by a whitespace.
   * `--dir`: the path of the _Quantification_ directory
   * `--channel_info_path`: the path to the .txt or .tsv file containing information on the image channels. The values must be separated by tabs. In particular, the file must contain a _Channel\_name_ and a _Cellular\_location_ (Nucleus or Cytoplasm) column.

   <details>
   <summary>Some other parameters may be required in case of non-standard uses of the pipeline, so they can usually be ignored</summary>

   * `--deepcell_path`: the absolute path to the Singularity image of DeepCell. The default is `/mnt/data2/shared/Lymphomoid-IF-software/deepcell.sif` where it's already present.
   * `--nucleus_channel`: the name of the channel associated to the nuclear marker (e.g. DAPI). The default is `DAPI`.

   </details>

An example of the script calling is: `python3 cellDetection.py --dir /mnt/data2/varrone/elisa_lymphomoids/HLS_Quantification/ --sample_names HLS25_7 HLS25_41acq01 HLS25_41acq03 --channel_info_path /mnt/data2/varrone/elisa_lymphomoids/mouse_channels.txt`

A lot of warning messages will appear, but they are normal. As long as the message `Channels extracted successfully.` appears, the software will have worked successfully. 

The run may take a while for each image (tens of minutes). For this reason, it is suggested to pass many images at the same time through the `--sample_names` parameters and run the script overnight.

### Marker quantification (i.e. segmentation)
This step obtains for each cell, from its mask detected in the cell detection step, the mean intensities of each of the markers in the nucleus and in the cytoplasm.

1. If the virtual environment has not been activated from the cell detection part, run the `source /mnt/data2/shared/Lymphomoid-IF-software/Lymphomoid-IF-venv/bin/activate` command to activate the virtual environment.

2. Then, run the _quantifyIntensities.py_ script, with the same required parameters as in the cell detection step:
   * `--sample_names`: the name of the images for cell detection, it could be 1 or many, separated by a whitespace.
   * `--dir`: the absolute path of the _HLS\_Quantification_ directory

   <details>
   <summary>Parameters for non-standard uses of the pipeline</summary>

   * `--nextflow_dir`: path to the Nextflow software. The default is `/mnt/data2/shared/Lymphomoid-IF-software/nextflow` where it's already present.

   </details>

An example of the script calling is: `python3 quantifyIntensities.py --dir /mnt/data2/varrone/elisa_lymphomoids/HLS_Quantification/ --sample_names HLS25_7 HLS25_41acq01 HLS25_41acq03`.


### Optional: downloading the software
As previously mentioned, the directory containing all the necessary software is already present in the uporicchiosrv1 server. 

If you want to download or update the four software required:
* bftools: bftools can be downloaded from the Bio-Formats [website](https://docs.openmicroscopy.org/bio-formats/6.8.1/users/comlinetools/index.html). Current version: 6.8.1.
* Virtual Environment: the Python packages required for running the pipeline are listed in the file _requirements.txt_. If you want to install the packages in your current virtual environment you can run `pip install -r requirements.txt`.
* DeepCell: the image for Singularity can be downloaded running `singularity pull deepcell.sif docker://vanvalenlab/deepcell-applications:latest`. Current version: 0.3.1.
* Nextflow: select the directory where ypu want to download nextflow and run `curl -s https://get.nextflow.io | bash`. For more information visit the [website](https://nextflow.io). Current version: 21.10.6.5660.
* MCMICRO: you can get the latest version of MCMICRO by moving to the directory where _nextflow_ is (e.g. `/mnt/data2/shared/Lymphomoid-IF-software/nextflow`) and running `./nextflow pull labsyspharm/mcmicro`. For more information visit the [website](https://mcmicro.org). Current version: Github revision 46abd97bc0.

## Step 5 - Classify cells and downstream analyses
In this step we will take quantification performed at step 4 and classify the cells according to the thresholds we calibrated at step 3. Each detected cell has a value for each of the cytoplasmic markers; if for a given cell they are all below the thresholds, that cell will be assigned to "otherCell", otherwise it will be assigned to the marker with the highest difference between the cytoplasmic intensity of that marker and its threshold. This means that in the current implementation a cell can be assigned to (or 'positive for') at most one marker. Additionally, cells will be classified as proliferating or not according to whether the Ki67 intensity is above or below its calibrated threshold. The script _ClassifyCells_Analyses.R_ performs this cell classification, downstream analyses and plotting for all lymphomoids provided as input. In detail:
1. Open _ClassifyCells_Analyses.R_ with a text editor or RStudio
2. Set up the input:
   * MainDir: your main directory, same as for the previous steps
   * ConfigTable: your configuration table, same as for the previous steps
   * Lymphomoids_to_process: set to "all" if you want all previously drawn lymphomoids to be processed, otherwise specify a subset of them as a character vector
   * Optionally, you can change some of the plotting parameters
3. If they are not already installed, install R packages `sp`, `ggplot2` and `reshape2`
4. Run the script

The script will generate the following output under MainDir (all folders are created automatically):
* `Classified_cells_tables/Table_<ImageName>_<PatientLymphomoidName>_AllCells.txt` containing, for each cell (rows), information on marker intensities, centroid coordinates in micrometers, assigned marker/antibody and wheter the cell is proliferating
* `Digital_IF_images/`
   - `IFimage_<ImageName>_<PatientLymphomoidName>_all.pdf` digitalized image with cells rendered as simple dots, color-coded by their assigned marker, and a red contour for the lymphomoid boundary
   - `IFimage_<ImageName>_<PatientLymphomoidName>_NoOtherCells_OnlyInLymphomoid.pdf` same but excluding otherCells and showing only cells inside the lymphomoid
   - `IFimage_<ImageName>_<PatientLymphomoidName>_NoOtherCells_OnlyInLymphomoid_OnlyProliferating.pdf` same but showing only proliferating cells inside the lymphomoid
* `*LymphomoidLevel*` files contain results summarized at the lymphomoid level, i.e. if two images contained the same lymphomoid (having the same name), their cell counts were averaged
   - `SummaryTable_LymphomoidLevel_AllCells.RData/.txt` table containing for each lymphomoid (rows), total cell counts for each marker, total counts of proliferating cells for each marker, and total number of cells detected in the lymphomoid
   - `LymphomoidLevel_CellTypeProportions_StackedBarplot.pdf/.txt` plot and table of cell type proportions across lymphomoids
   - `LymphomoidLevel_CellTypeProportions_ExclOtherCells_StackedBarplot.pdf/.txt` same as above but excluding otherCells
   - `LymphomoidLevel_ProlifVsNot_<marker>.pdf/.txt` plot and table of proportion of proliferating/not proliferating cells of a given marker across lymphomoids
* `Results_ImageXLymphomoid_level/*` same output of the point above, but without lymphomoid-level summarization (i.e. 'image X lymphomoid' level). It could be useful to check that two images of the same lymphomoid (before summarization) have comparable cell proportions
* `log_files/*` .log files with date and time of the script runs in their name and containing the various input file paths for reproducibility

Since all tables are saved, you can customize your plots/perfom additional downstream analyses by reading them back in your favourite programming language.



