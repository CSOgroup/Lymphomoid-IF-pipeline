
println "Script starting..."

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.RectangleROI
import qupath.lib.gui.scripting.QPEx

def imageData = QPEx.getCurrentImageData()
clearDetections();
selectAnnotations();

// ----------- Input ----------- 
OutDir = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Processed/Pipeline_test/" // Absolute path to your desired output directory
ImageName = "HLS33_21acq03"
DAPI_thresh = 200
FITC_thresh = 300
CY5_thresh = 100
CFP_thresh = 300
RFP_thresh = 300
Alexa_thresh = 300
// -----------------------------


// ----------- Calibration ----------- 
// DAPI
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 15.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.0,  "minAreaMicrons": 5.0,  "maxAreaMicrons": 100.0,  "threshold": '+DAPI_thresh+',  "watershedPostProcess": true,  "cellExpansionMicrons": 2.5,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
// // FITC (B-cells)
// setCellIntensityClassifications("Cytoplasm: FITC mean", FITC_thresh)
// // CY5 (Ki67)
// setCellIntensityClassifications("Nucleus: CY5 mean", CY5_thresh)
// // CFP (CD4 T-cells)
// setCellIntensityClassifications("Cytoplasm: CFP mean", CFP_thresh)
// // RFP (CD8 T-cells)
// setCellIntensityClassifications("Cytoplasm: RFP mean", RFP_thresh)
// // Alexa 594 (macrophages)
// setCellIntensityClassifications("Cytoplasm: Alexa 594 mean", Alexa_thresh)
// -----------------------------------


// ----------- Saving ----------- 
File file = new File(OutDir+ImageName+"_AllThresholds.txt")
file.write "DAPI_thresh "+DAPI_thresh+"\n"+"FITC_thresh "+FITC_thresh+"\n"+"CY5_thresh "+CY5_thresh+"\n"+"CFP_thresh "+CFP_thresh+"\n"+"RFP_thresh "+RFP_thresh+"\n"+"Alexa_thresh "+Alexa_thresh+"\n"
// ------------------------------ 

println "Done"
