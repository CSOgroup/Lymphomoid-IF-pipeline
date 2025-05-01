
println "Script starting..."

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.RectangleROI
import qupath.lib.gui.scripting.QPEx

def imageData = QPEx.getCurrentImageData()

// ----------- Input ----------- 
MainDir = "/Users/daniele/Mounted_folder8/Daniele/elisa_lymphomoids/test2/" // Absolute path to your main directory
ConfigTable = "/Users/daniele/Mounted_folder8/Daniele/elisa_lymphomoids/Lymphomoid-IF-pipeline/human_channels.txt" // Configuration table as 'mouse_channels.txt' or 'human_channels.txt' 
ImageName = "expand_to_all_available_images" // Either a specific 'ImageName' (e.g. 'HLS16_s01_acq01'), or 'expand_to_all_available_images'. In the latter case, the same calibration table will be copied for all images under the 'Lymphomoid_boundaries/' folder (careful, not calibrating the thresholds separately for each image might result in unrealistic cell type calls).
DAPI_thresh = 100
FITC_thresh = 280
CY5_thresh = 110
CFP_thresh = 280
RFP_thresh = 300
Alexa_thresh = 100
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
// // DAPI
// print('DAPI: detecting '+channel_antibody_map['DAPI']+' in the '+channel_celllocation_map['DAPI'])
// setCellIntensityClassifications("Cell: DAPI: Mean", DAPI_thresh)
// // FITC
// print('FITC: detecting '+channel_antibody_map['FITC']+' in the '+channel_celllocation_map['FITC'])
// setCellIntensityClassifications("Cell: FITC: Mean", FITC_thresh)
// // CY5
// print('CY5: detecting '+channel_antibody_map['CY5']+' in the '+channel_celllocation_map['CY5'])
// setCellIntensityClassifications("Cell: CY5: Mean", CY5_thresh)
// // CFP
// print('CFP: detecting '+channel_antibody_map['CFP']+' in the '+channel_celllocation_map['CFP'])
// setCellIntensityClassifications("Cell: CFP: Mean", CFP_thresh)
// // RFP
// print('RFP: detecting '+channel_antibody_map['RFP']+' in the '+channel_celllocation_map['RFP'])
// setCellIntensityClassifications("Cell: RFP: Mean", RFP_thresh)
// // Alexa 594
// print('Alexa 594: detecting '+channel_antibody_map['Alexa 594']+' in the '+channel_celllocation_map['Alexa 594'])
// setCellIntensityClassifications("Cell: Alexa 594: Mean", Alexa_thresh)
// -----------------------------------


// ----------- Saving ----------- 
new File(MainDir+"/Calibrated_thresholds/").mkdirs()
if (ImageName=="expand_to_all_available_images"){
    dh = new File(MainDir+'/Lymphomoid_boundaries/')
    dh.eachFile {
        String thiss = it
        thiss = thiss.replaceAll(".+/", "");
        thiss = thiss.substring(0,thiss.lastIndexOf("_"));
        thiss = thiss.substring(0,thiss.lastIndexOf("_"));
        File file = new File(MainDir+"/Calibrated_thresholds/"+thiss+"_AllThresholds.txt")
        file.write "DAPI_thresh "+DAPI_thresh+"\n"+"FITC_thresh "+FITC_thresh+"\n"+"CY5_thresh "+CY5_thresh+"\n"+"CFP_thresh "+CFP_thresh+"\n"+"RFP_thresh "+RFP_thresh+"\n"+"Alexa_thresh "+Alexa_thresh+"\n"        
    }
} else {
    File file = new File(MainDir+"/Calibrated_thresholds/"+ImageName+"_AllThresholds.txt")
    file.write "DAPI_thresh "+DAPI_thresh+"\n"+"FITC_thresh "+FITC_thresh+"\n"+"CY5_thresh "+CY5_thresh+"\n"+"CFP_thresh "+CFP_thresh+"\n"+"RFP_thresh "+RFP_thresh+"\n"+"Alexa_thresh "+Alexa_thresh+"\n"        
}
println "Calibrated thresholds saved"
// ------------------------------ 

println "Done"