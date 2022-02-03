
println "Script starting..."

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.gui.scripting.QPEx

def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
double pixelDimMicrons = server.getPixelCalibration().getAveragedPixelSizeMicrons()

// ----------- Input ----------- 
MainDir = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Processed/Pipeline_final_test/" // Absolute path to your main directory
ImageName = "HLS16_s01_acq02"
PatientLymphomoidName = "LP02ctrl01" // Avoid underscores for PatientLymphomoidName
// -----------------------------

// ------ Saving boundary ------ 
OutFile = ImageName+"_"+PatientLymphomoidName+"_Boundary.txt" // Output file name 
if (!getSelectedROI()) { println "No boundary selected" }
def polyPoints = getSelectedROI().getAllPoints()
String result = polyPoints.join(";")
new File(MainDir+"Lymphomoid_boundaries/").mkdirs()
new File(MainDir+"Tiff/").mkdirs()
new File(MainDir+"Quantification/"+ImageName+"/registration").mkdirs()
File file = new File(MainDir+"Lymphomoid_boundaries/"+OutFile)
file.write result+"\n"
file << "PixelDim="+pixelDimMicrons.toString()+"\n"
// -----------------------------
println "Boundary of "+ImageName+", "+PatientLymphomoidName+" saved."
println "Done"
