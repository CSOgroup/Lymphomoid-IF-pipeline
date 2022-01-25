
println "Script starting..."

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.gui.scripting.QPEx

def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
double pixelDimMicrons = server.getPixelCalibration().getAveragedPixelSizeMicrons()

// ----------- Input ----------- 
OutDir = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Processed/Pipeline_test/" // Absolute path to your desired output directory
ImageName = "HLS33_s21_acq03"
LymphomoidName = "Pembroluzimab01" // Avoid underscores for LymphomoidName
// -----------------------------

OutFile = ImageName+"_"+LymphomoidName+"_Boundary.txt" // Output file name
def polyPoints = getSelectedROI().getAllPoints()
String result = polyPoints.join(";")
File file = new File(OutDir+OutFile)
file.write result+"\n"
file << "PixelDim="+pixelDimMicrons.toString()+"\n"

println "Done"
