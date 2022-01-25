
println "Script starting..."

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.gui.scripting.QPEx

def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
double pixelDimMicrons = server.getPixelCalibration().getAveragedPixelSizeMicrons()

// ----------- Input ----------- 
MainDir = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Processed/Pipeline_test/" // Absolute path to your main directory
ImageName = "HLS01_s02_acq03"
LymphomoidName = "Pembroluzimab02" // Avoid underscores for LymphomoidName
// -----------------------------

// ------ Saving boundary ------ 
OutFile = ImageName+"_"+LymphomoidName+"_Boundary.txt" // Output file name 
def polyPoints = getSelectedROI().getAllPoints()
String result = polyPoints.join(";")
new File(MainDir+"Lymphomoid_boundaries/").mkdirs()
File file = new File(MainDir+"Lymphomoid_boundaries/"+OutFile)
file.write result+"\n"
file << "PixelDim="+pixelDimMicrons.toString()+"\n"
// -----------------------------

println "Done"
