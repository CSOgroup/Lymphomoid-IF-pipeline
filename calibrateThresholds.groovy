
println "Script starting..."

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.RectangleROI
import qupath.lib.gui.scripting.QPEx

def imageData = QPEx.getCurrentImageData()
clearDetections();

// ----------- Input ----------- 
MainDir = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Processed/Pipeline_test/" // Absolute path to your main directory
ConfigTable = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Scripts/Lymphomoid-IF-pipeline/human_channels.txt" // Configuration table as 'mouse_channels.txt' or 'human_channels.txt' 
ImageName = "HLS01_s02_acq03"
DAPI_thresh = 200
FITC_thresh = 300
CY5_thresh = 110
CFP_thresh = 300
RFP_thresh = 300
Alexa_thresh = 300
// -----------------------------

// ------ Parsing config table ------
def fileArray = new File(ConfigTable) as String[]
channel_celllocation_map = [:]
channel_antibody_map = [:]
for (int line in 1..(fileArray.size()-1)){ 
    channel_celllocation_map[fileArray[line].split("\t", -1)[0]] = fileArray[line].split("\t", -1)[1]
    channel_antibody_map[fileArray[line].split("\t", -1)[0]] = fileArray[line].split("\t", -1)[2]
}
// ----------------------------------


// ----------- Calibration ----------- 
// DAPI
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 15.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.0,  "minAreaMicrons": 5.0,  "maxAreaMicrons": 100.0,  "threshold": '+DAPI_thresh+',  "watershedPostProcess": true,  "cellExpansionMicrons": 2.5,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
// // FITC
// print('FITC: detecting '+channel_antibody_map['FITC']+' in the '+channel_celllocation_map['FITC'])
// setCellIntensityClassifications(channel_celllocation_map['FITC']+": FITC mean", FITC_thresh)
// // CY5
// print('CY5: detecting '+channel_antibody_map['CY5']+' in the '+channel_celllocation_map['CY5'])
// setCellIntensityClassifications(channel_celllocation_map['CY5']+": CY5 mean", CY5_thresh)
// // CFP
// print('CFP: detecting '+channel_antibody_map['CFP']+' in the '+channel_celllocation_map['CFP'])
// setCellIntensityClassifications(channel_celllocation_map['CFP']+": CFP mean", CFP_thresh)
// // RFP
// print('RFP: detecting '+channel_antibody_map['RFP']+' in the '+channel_celllocation_map['RFP'])
// setCellIntensityClassifications(channel_celllocation_map['RFP']+": RFP mean", RFP_thresh)
// // Alexa 594
// print('Alexa 594: detecting '+channel_antibody_map['Alexa 594']+' in the '+channel_celllocation_map['Alexa 594'])
// setCellIntensityClassifications(channel_celllocation_map['Alexa 594']+": Alexa 594 mean", Alexa_thresh)
// -----------------------------------


// ----------- Saving ----------- 
new File(MainDir+"Calibrated_thresholds/").mkdirs()
File file = new File(MainDir+"Calibrated_thresholds/"+ImageName+"_AllThresholds.txt")
file.write "DAPI_thresh "+DAPI_thresh+"\n"+"FITC_thresh "+FITC_thresh+"\n"+"CY5_thresh "+CY5_thresh+"\n"+"CFP_thresh "+CFP_thresh+"\n"+"RFP_thresh "+RFP_thresh+"\n"+"Alexa_thresh "+Alexa_thresh+"\n"
// ------------------------------ 

println "Done"
