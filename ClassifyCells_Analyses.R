
########## Input ##########
MainDir = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Processed/Pipeline_test/" # Absolute path to your main directory
ConfigTable = "/Users/daniele/Mounted_folder_daniele_ndata/elisa_lymphomoids/Scripts/Lymphomoid-IF-pipeline/mouse_channels.txt" # Configuration table as 'mouse_channels.txt' or 'human_channels.txt' 
Lymphomoids_to_process = "all" # "all", or vector of boundary file names (e.g. Lymphomoids_to_process = c( "HLS01_s02_acq03_Pembroluzimab01_Boundary.txt", "HLS01_s02_acq03_Pembroluzimab02_Boundary.txt" ) )
###########################

###### Plotting parameters ######
plot_IF_images = TRUE # TRUE or FALSE. FALSE makes the script run faster.
antibody_colors = data.frame(antibody_mouse = c("B220","CD4","CD8","F4/80","otherCell"), antibody_human = c("CD20","CD4","CD8","CD68","otherCell"), color = c("green","yellow","cyan","magenta","gray"))
maxSizePlot_inches = 20
cellSizePlot_um = 100
boundary_thickness = 2
prolif_color = "deepskyblue"
#################################

########## Functions ##########
library(sp)
library(ggplot2)
library(reshape2)

parse_boundary = function(pl, pixelSize){
   pl_points = unlist(strsplit(as.character(pl[1,"V1"]), split=";"))
   plp_df = data.frame(x = rep(NA,length(pl_points)+1), y = rep(NA,length(pl_points)+1))
   index = 1
   for (point in pl_points)
   {
      poiSpl = unlist(strsplit(point, split = ", "))
      plp_df[index,"x"] = as.numeric(substr(poiSpl[1],8,nchar(poiSpl[1])))*pixelSize
      plp_df[index,"y"] = (as.numeric((poiSpl[2])))*pixelSize
      index = index+1
   }
   plp_df[nrow(plp_df),"x"] = plp_df[1,"x"]
   plp_df[nrow(plp_df),"y"] = plp_df[1,"y"]
   return(plp_df)
}

classify_cells = function(quant, ll_config){
   this_inte = sweep(quant[,rownames(ll_config)[!(rownames(ll_config) %in% c( "DAPI","Ki67" ))] ],2,ll_config[rownames(ll_config)[!(rownames(ll_config) %in% c( "DAPI","Ki67" ))],"calibrated_thresholds"])
   is_cell = rowSums(this_inte>0)>0
   this_inte = this_inte[is_cell %in% c(T),]
   inteMax = apply(this_inte,1,which.max)
   quant$spatial_1 = quant$X_centroid_nucleus*pixelSize
   quant$spatial_2 = quant$Y_centroid_nucleus*pixelSize
   quant$CellType_index = 0
   quant[names(inteMax),"CellType_index"] = as.character(inteMax)
   quant$CellType_color = "gray"
   quant$CellType_antibody = "otherCell"
   quant$CellType_marker = "otherCell"
   for (index in 1:4)
   {
      quant[quant$CellType_index==index,"CellType_color"] = antibody_colors[colnames(this_inte)[index],"color"]
      quant[quant$CellType_index==index,"CellType_antibody"] = colnames(this_inte)[index]
      quant[quant$CellType_index==index,"CellType_marker"] = ll_config[colnames(this_inte)[index],"Marker_of"]
   }
   quant$CellType_marker = factor(quant$CellType_marker, levels = c("otherCell",ll_config[colnames(this_inte),"Marker_of"]))
   quant$CellType_antibody = factor(quant$CellType_antibody, levels = c("otherCell",colnames(this_inte)))
   quant$CellType_index = NULL
   quant[,"is_proliferating"] = as.numeric(quant[,"Ki67"])>ll_config["Ki67","calibrated_thresholds"]
   return(quant)
}

plot_digital_image = function(fileName, quant, plp_df, withOtherCells = T, onlyInLymphomoid = F, onlyProliferating = F){
   if (onlyProliferating) { quant = quant[quant$is_proliferating,] }
   if (!withOtherCells) { quant = quant[quant$CellType_antibody != "otherCell",] }
   if (onlyInLymphomoid) { quant = quant[quant$in_lymphomoid,] }
   abdf = data.frame(row.names = as.character(unique(quant$CellType_antibody)), ab = as.character(unique(quant$CellType_antibody)), color = unique(quant$CellType_color))
   colorz = abdf[intersect(levels(quant$CellType_antibody),rownames(abdf)),"color"]
   widthz = max(c(quant$spatial_1,plp_df$x))-min(c(quant$spatial_1,plp_df$x))
   heightz = max(c(quant$spatial_2,plp_df$y))-min(c(quant$spatial_2,plp_df$y))
   rescaling_to_plot = maxSizePlot_inches/max(c(widthz,heightz))
   r = cellSizePlot_um * rescaling_to_plot
   pdf( file = fileName, width = widthz*rescaling_to_plot, height = heightz*rescaling_to_plot, useDingbats = F)
   p = ggplot(quant, aes(x=spatial_1, y=spatial_2, color=CellType_antibody)) + geom_point(stroke=0,size=r) + scale_color_manual(values=colorz ) + scale_y_reverse() + geom_path(data = plp_df, mapping = aes(x = x, y = y), color = "red", size = boundary_thickness) + theme_void() + theme(panel.background = element_rect(fill = 'black', colour = 'black')) + theme(legend.position = 'none')
   print(p)
   dev.off()
   return(p)
}

plot_CellTypeProportions = function(markers, ldf, OutFileRoot){
   colors = antibody_colors[markers,"color"]
   cdf = ldf[,markers]/apply(ldf[,markers],1,sum)
   cdf$Sample = rownames(cdf)
   write.table(cdf, file = paste0(MainDir,OutFileRoot,".txt"), quote = F, col.names = T, row.names = T, sep = "\t")
   cdf = melt(cdf,id="Sample")
   cdf$variable = factor(cdf$variable, levels = markers)
   cdf$Sample = factor(cdf$Sample, levels = rownames(ldf))
   pdf(paste0(MainDir,OutFileRoot,".pdf"),2*nrow(ldf),6)
   p = ggplot(data=cdf, aes(x=Sample, y=value, fill=variable)) +
      geom_bar(stat="identity", colour="black", size = 0.5) + ylab("Proportion") + xlab("") + scale_fill_manual(name = "Marker",values=colors) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(labels= rownames(ldf))
   print(p)
   dev.off()
   return(p)
}

plot_ProlifVsNot = function(markers, ldf, OutFileRoot){
   colors = antibody_colors[markers,"color"]
   for (marker in markers)
   {
      cdf = ldf[,c(marker,paste0("prolif_",marker))]
      total = cdf[,marker]
      cdf[,marker] = cdf[,marker]-cdf[,paste0("prolif_",marker)]
      cdf = cdf/total
      colnames(cdf) = c("Not proliferating","Proliferating")
      cdf$Sample = rownames(cdf)
      write.table(cdf, file = paste0(MainDir,OutFileRoot,gsub("/","-",marker),".txt"), quote = F, col.names = T, row.names = T, sep = "\t")
      cdf = melt(cdf,id="Sample")
      cdf$variable = factor(cdf$variable, levels = c("Not proliferating","Proliferating"))
      cdf$Sample = factor(cdf$Sample, levels = rownames(ldf))
      pdf(paste0(MainDir,OutFileRoot,gsub("/","-",marker),".pdf"),2*nrow(ldf),6)
      p = ggplot(data=cdf, aes(x=Sample, y=value, fill=variable)) +
         geom_bar(stat="identity", colour="black", size = 0.5) + ylab("Proportion") + xlab("") + scale_fill_manual(name = marker,values=c(colors[which(marker==markers)],prolif_color)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(labels= rownames(ldf))
      print(p)
      dev.off()
   }
   return(p)
}

###############################

########## Main ##########
channels = read.table( file = ConfigTable, sep = "\t", header = T, quote = '' ,stringsAsFactors = F)
dir.create(paste0(MainDir,"Digital_IF_images/"), showWarnings = F)
dir.create(paste0(MainDir,"Classified_cells_tables/"), showWarnings = F)
if (Lymphomoids_to_process=="all"){ Lymphomoids_to_process = list.files(paste0( MainDir,"Lymphomoid_boundaries/" ), pattern = "Boundary.txt") }
if (length(intersect(channels$Antibody,antibody_colors$antibody_mouse))==4 ) rownames(antibody_colors) = antibody_colors$antibody_mouse
if (length(intersect(channels$Antibody,antibody_colors$antibody_human))==4 ) rownames(antibody_colors) = antibody_colors$antibody_human
all_markerz = c("otherCell", channels$Antibody[!(channels$Antibody %in% c( "DAPI","Ki67" ))])
all_colz = c("ImageName","PatientLymphomoidName",all_markerz, paste0("prolif_",all_markerz ),"TotalCells")
tdf = data.frame(matrix(nrow = 0, ncol = length(all_colz), dimnames = list(NULL,all_colz)))
colnames(tdf)[colnames(tdf)=="F4.80"] = "F4/80"
colnames(tdf)[colnames(tdf)=="prolif_F4.80"] = "prolif_F4/80"
quant_file_vec = c() # for log file
calib_file_vec = c() # for log file
for (ll in Lymphomoids_to_process)
{
   cat("\n","Processing",ll,"...","\n" )
   ## Parsing and reading
   ImageName = sub("_[^_]+$", "", substr(ll,1,nchar(ll)-13))
   PatientLymphomoidName = substr(ll,nchar(ImageName)+2,nchar(ll)-13)
   quant_file = paste0(MainDir,"HLS_Quantification/",ImageName,"/quantification/mesmer-",ImageName,"_merged.csv")
   quant_file_vec = c(quant_file_vec,substr(quant_file,nchar(MainDir)+1,nchar(quant_file)))
   if (!file.exists(quant_file)) { next }
   quant = read.csv(quant_file,stringsAsFactors = F)
   quant = quant[,!(colnames(quant) %in% c( "Area_nucleus","MajorAxisLength_nucleus","MinorAxisLength_nucleus","Eccentricity_nucleus","Solidity_nucleus","Extent_nucleus","Orientation_nucleus","X_centroid_cytoplasm","Y_centroid_cytoplasm","Area_cytoplasm","MajorAxisLength_cytoplasm","MinorAxisLength_cytoplasm","Eccentricity_cytoplasm","Solidity_cytoplasm","Extent_cytoplasm","Orientation_cytoplasm" ))]
   colnames(quant)[colnames(quant)=="F4.80"] = "F4/80"
   calib_file = paste0(MainDir,"Calibrated_thresholds/",ImageName,"_AllThresholds.txt")
   calib_file_vec = c(calib_file_vec,substr(calib_file,nchar(MainDir)+1,nchar(calib_file)))
   if (!file.exists(calib_file)) { next }
   calib = read.table( file = calib_file, sep = " ", header = F, quote = '' ,stringsAsFactors = F)
   rownames(calib) = calib$V1
   pl = read.table(file = paste0(MainDir,"Lymphomoid_boundaries/",ll), sep = "-")
   pixelSize = as.numeric(substr(as.character(pl[2,"V1"]),10,nchar(as.character(pl[2,"V1"]))))
   plp_df = parse_boundary(pl,pixelSize)

   ## Setting up config table of this lymphomoid
   ll_config = channels
   rownames(ll_config) = paste0(sub(" .*", "", ll_config$Channel_name),"_thresh")
   ll_config$calibrated_thresholds = calib[rownames(ll_config),"V2"]
   rownames(ll_config) = ll_config$Antibody

   ## Classifying cells and checking which ones are inside the lymphomoid
   quant = classify_cells(quant, ll_config)
   inPolygon = point.in.polygon(point.x = quant$spatial_1, point.y = quant$spatial_2, pol.x = plp_df$x, pol.y = plp_df$y)
   quant$in_lymphomoid = inPolygon==1

   ## Plotting digital IF images
   if (plot_IF_images)
   {
      fileName = paste0(MainDir,"Digital_IF_images/","IFimage_",ImageName,"_",PatientLymphomoidName,"_all.pdf")
      p = plot_digital_image(fileName, quant, plp_df, withOtherCells = T, onlyInLymphomoid = F, onlyProliferating = F)
      fileName = paste0(MainDir,"Digital_IF_images/","IFimage_",ImageName,"_",PatientLymphomoidName,"_NoOtherCells_OnlyInLymphomoid.pdf")
      p = plot_digital_image(fileName, quant, plp_df, withOtherCells = F, onlyInLymphomoid = T, onlyProliferating = F)
      fileName = paste0(MainDir,"Digital_IF_images/","IFimage_",ImageName,"_",PatientLymphomoidName,"_NoOtherCells_OnlyInLymphomoid_OnlyProliferating.pdf")
      p = plot_digital_image(fileName, quant, plp_df, withOtherCells = F, onlyInLymphomoid = T, onlyProliferating = T)
   }
   
   ## Saving data
   quant$X_centroid_nucleus = NULL
   quant$Y_centroid_nucleus = NULL
   colnames(quant)[colnames(quant) %in% c( "spatial_1","spatial_2" )] = c( "x_centroid_nucleus_um","y_centroid_nucleus_um" )
   quant = quant[quant$in_lymphomoid,]
   write.table(quant, file = paste0(MainDir,"Classified_cells_tables/Table_",ImageName,"_",PatientLymphomoidName,"_AllCells.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

   ## Concatenating summary for all lymphomoids
   ttdf = data.frame(matrix(0,nrow = 1, ncol = length(all_colz), dimnames = list(paste0(ImageName,"_",PatientLymphomoidName),all_colz)))
   colnames(ttdf)[colnames(ttdf)=="F4.80"] = "F4/80"
   colnames(ttdf)[colnames(ttdf)=="prolif_F4.80"] = "prolif_F4/80"
   ttdf[, c( "ImageName","PatientLymphomoidName" )] = c( ImageName,PatientLymphomoidName )
   ttdf[,names(table(quant$CellType_antibody))] = as.numeric(table(quant$CellType_antibody))
   ttdf[,"TotalCells"] = sum(as.numeric(table(quant$CellType_antibody)))
   ttdf[,paste0("prolif_",rownames(table(quant$CellType_antibody,quant$is_proliferating)))] = as.numeric(table(quant$CellType_antibody,quant$is_proliferating)[,"TRUE"])
   tdf = rbind(tdf, ttdf)
   
}

## Saving and plotting at the image X lymphomoid level
dir.create(paste0(MainDir,"Results_ImageXLymphomoid_level/"),showWarnings = F)
save(tdf, file = paste0(MainDir,"Results_ImageXLymphomoid_level/SummaryTable_ImageXLymphomoid_AllCells.RData"))
write.table(tdf, file = paste0(MainDir,"Results_ImageXLymphomoid_level/SummaryTable_ImageXLymphomoid_AllCells.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
p = plot_CellTypeProportions(markers = all_markerz, ldf = tdf, OutFileRoot = "Results_ImageXLymphomoid_level/ImageXLymphomoidLevel_CellTypeProportions_StackedBarplot")
p = plot_CellTypeProportions(markers = all_markerz[all_markerz!="otherCell"], ldf = tdf, OutFileRoot = "Results_ImageXLymphomoid_level/ImageXLymphomoidLevel_CellTypeProportions_ExclOtherCells_StackedBarplot")
p = plot_ProlifVsNot(markers = all_markerz[all_markerz!="otherCell"], ldf = tdf, OutFileRoot = "Results_ImageXLymphomoid_level/ImageXLymphomoidLevel_ProlifVsNot_")

## Consolidating the same lymphomoid across images, saving and plotting
tdf$ImageName = NULL
ldf = aggregate(.~PatientLymphomoidName, data = tdf, FUN=mean)
rownames(ldf) = ldf$PatientLymphomoidName
save(ldf, file = paste0(MainDir,"SummaryTable_LymphomoidLevel_AllCells.RData"))
write.table(ldf, file = paste0(MainDir,"SummaryTable_LymphomoidLevel_AllCells.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
p = plot_CellTypeProportions(markers = all_markerz, ldf = ldf, OutFileRoot = "LymphomoidLevel_CellTypeProportions_StackedBarplot")
p = plot_CellTypeProportions(markers = all_markerz[all_markerz!="otherCell"], ldf = ldf, OutFileRoot = "LymphomoidLevel_CellTypeProportions_ExclOtherCells_StackedBarplot")
p = plot_ProlifVsNot(markers = all_markerz[all_markerz!="otherCell"], ldf = ldf, OutFileRoot = "LymphomoidLevel_ProlifVsNot_")

## Creating and saving log table
logdf = data.frame(
   MainDir = MainDir,
   ConfigTable = ConfigTable,
   lymphomoid_boundaries = paste(paste0("Lymphomoid_boundaries/",Lymphomoids_to_process), collapse = ","),
   quantification_tables = paste(quant_file_vec, collapse = ","),
   calibration_tables = paste(calib_file_vec, collapse = ","),
   ImageXLymphomoid_level_output_table = paste0("Results_ImageXLymphomoid_level/SummaryTable_ImageXLymphomoid_AllCells.txt"),
   Lymphomoid_level_output_table = paste0("SummaryTable_LymphomoidLevel_AllCells.txt")
   )
tl = data.frame(t(logdf))
dir.create(paste0(MainDir,"log_files/"),showWarnings=F)
logdf = data.frame(variable = paste0(rownames(tl)," = "), value = tl[,1])
logfilename = paste0("run_ClassifyCells_date",Sys.Date(),"_time",gsub(":",".",substr(Sys.time(),nchar(as.character(Sys.Date()))+2,nchar(as.character(Sys.time())))),".log")
write.table(logdf, file = paste0(MainDir,"log_files/",logfilename), row.names = F, col.names = F, sep = "\t", quote = F)


