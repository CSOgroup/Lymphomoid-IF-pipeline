
########## Input ##########
MainDir = "/mnt/ndata/daniele/elisa_lymphomoids/Processed/Pipeline_test/" # Absolute path to your main directory
ConfigTable = "/mnt/ndata/daniele/elisa_lymphomoids/Lymphomoid-IF-pipeline/mouse_channels.txt" # Configuration table as 'mouse_channels.txt' or 'human_channels.txt' 
Lymphomoids_to_process = "all" # "all", or vector of boundary file names (e.g. Lymphomoids_to_process = c( "HLS01_s02_acq03_Pembroluzimab01_Boundary.txt", "HLS01_s02_acq03_Pembroluzimab02_Boundary.txt" ) )
###########################

###### Plotting parameters ######
plot_IF_images = TRUE # TRUE or FALSE. FALSE makes the script run faster.

# ToDo: there could be arbitrary channels, so use a palette
antibody_colors = data.frame(antibody_mouse = c("B220","CD4","CD8","F4/80","otherCell"), antibody_human = c("CD20","CD4","CD8","CD68","otherCell"), color = c("green","yellow","cyan","magenta","gray"), stringsAsFactors = F)
maxSizePlot_inches = 20
cellSizePlot_um = 100
boundary_thickness = 2
#################################

########## Functions ##########
library(sp)
library(ggplot2)
library(reshape2)

parse_boundary = function(pl, pixelSize){
   # Parse the boundary file to get the coordinates of the boundary points
   pl_points = unlist(strsplit(as.character(pl[1,"V1"]), split=";"))

   # Initialize a data frame to store the coordinates
   plp_df = data.frame(x = rep(NA,length(pl_points)+1), y = rep(NA,length(pl_points)+1))
   index = 1

   # Loop through each point and extract the coordinates
   for (point in pl_points)
   {
      poiSpl = unlist(strsplit(point, split = ", "))
      plp_df[index,"x"] = as.numeric(substr(poiSpl[1],8,nchar(poiSpl[1])))*pixelSize
      plp_df[index,"y"] = (as.numeric((poiSpl[2])))*pixelSize
      index = index+1
   }

   # Close the loop by connecting the last point to the first point
   plp_df[nrow(plp_df),"x"] = plp_df[1,"x"]
   plp_df[nrow(plp_df),"y"] = plp_df[1,"y"]
   return(plp_df)
}

classify_cells = function(quant, ll_config){
   # Classify cells based on the thresholded intensity values only for the cell type markers
   this_inte = sweep(quant[,rownames(ll_config)[ll_config$Is_CellType_Marker]],2,ll_config[rownames(ll_config)[ll_config$Is_CellType_Marker],"calibrated_thresholds"])
   
   # Identify cells that have positive intensity values for any marker
   is_cell = rowSums(this_inte>0)>0

   # Filter the intensity data to only include cells that are classified as cells
   this_inte = this_inte[is_cell %in% c(T),]
   
   # Find the index of the maximum intensity value for each cell
   inteMax = apply(this_inte,1,which.max)
   
   # Add spatial coordinates and cell type indices to the quantification data
   quant$spatial_1 = quant$X_centroid_nucleus*pixelSize
   quant$spatial_2 = quant$Y_centroid_nucleus*pixelSize
   quant$CellType_index = 0
   quant[names(inteMax),"CellType_index"] = as.character(inteMax)
   quant$CellType_color = "gray"
   quant$CellType_antibody = "otherCell"
   quant$CellType_marker = "otherCell"
   
   # Assign colors, antibodies, and markers based on the maximum intensity index
   for (index in 1:4)
   {
      quant[quant$CellType_index==index,"CellType_color"] = antibody_colors[colnames(this_inte)[index],"color"]
      quant[quant$CellType_index==index,"CellType_antibody"] = colnames(this_inte)[index]
      quant[quant$CellType_index==index,"CellType_marker"] = ll_config[colnames(this_inte)[index],"Marker_of"]
   }
   
   quant$CellType_marker = factor(quant$CellType_marker, levels = c("otherCell",ll_config[colnames(this_inte),"Marker_of"]))
   quant$CellType_antibody = factor(quant$CellType_antibody, levels = c("otherCell",colnames(this_inte)))
   quant$CellType_index = NULL

   # For every non-cell type marker taken from the configuration table, add a column indicating whether the cell is positive or negative for that marker
   for (marker in ll_config[ll_config$Is_CellType_Marker == FALSE & ll_config$Antibody != "DAPI", "Antibody"])
   {
      quant[,paste0(marker,"+")] = quant[,marker]>ll_config[marker,"calibrated_thresholds"]
   }

   return(quant)
}

plot_digital_image = function(fileName, quant, plp_df, withOtherCells = T, onlyInLymphomoid = F, nonCellTypeMarker = NULL){
   if (!withOtherCells) { quant = quant[quant$CellType_antibody != "otherCell",] }
   if (onlyInLymphomoid) { quant = quant[quant$in_lymphomoid,] }
   # If a non-cell type marker is provided, plot only cells that have column "marker_positive" equal to TRUE
   if (!is.null(nonCellTypeMarker)) { quant = quant[quant[,paste0(nonCellTypeMarker,"+")]==TRUE,] }
   
   abdf = data.frame(row.names = as.character(unique(quant$CellType_antibody)), ab = as.character(unique(quant$CellType_antibody)), color = unique(quant$CellType_color), stringsAsFactors = F)
   colorz = abdf[intersect(levels(quant$CellType_antibody),rownames(abdf)),"color"]
   widthz = max(c(quant$spatial_1,plp_df$x))-min(c(quant$spatial_1,plp_df$x))
   heightz = max(c(quant$spatial_2,plp_df$y))-min(c(quant$spatial_2,plp_df$y))
   rescaling_to_plot = maxSizePlot_inches/max(c(widthz,heightz))
   r = cellSizePlot_um * rescaling_to_plot
   pdf( file = fileName, width = widthz*rescaling_to_plot, height = heightz*rescaling_to_plot, useDingbats = F)
   p = ggplot(quant, aes(x=spatial_1, y=spatial_2, color=CellType_antibody)) + geom_point(stroke=0,size=r) + scale_color_manual(values=colorz ) + scale_y_reverse() + geom_path(data = plp_df, mapping = aes(x = x, y = y), color = "red", size = boundary_thickness) + theme_void() + theme(legend.position = 'none') + theme(panel.background = element_rect(fill = 'black', colour = 'black'), plot.background = element_rect(fill = "black"))
   print(p)
   dev.off()
   return(p)
}

# Function to plot proportions of different cell types in each sample
plot_CellTypeProportions = function(markers, ldf, OutFileRoot){
   # Get colors for each marker from antibody_colors table
   colors = antibody_colors[markers,"color"]
   
   # Calculate proportions of each marker relative to total cells
   cdf = ldf[,markers]/apply(ldf[,markers],1,sum)
   cdf$Sample = rownames(cdf)
   
   # Write proportions to output file
   write.table(cdf, file = paste0(MainDir,OutFileRoot,".txt"), quote = F, col.names = T, row.names = T, sep = "\t")
   
   # Reshape data for plotting
   cdf = melt(cdf,id="Sample")
   cdf = merge(cdf,ldf[,c("PatientLymphomoidName","TotalCells")],by.x = "Sample", by.y = "row.names")
   cdf$variable = factor(cdf$variable, levels = markers)
   cdf$Sample = factor(cdf$Sample, levels = rownames(ldf))
   cdf$TotalCells = paste0("N=",round(cdf$TotalCells))
   
   # Create stacked bar plot
   pdf(paste0(MainDir,OutFileRoot,".pdf"),2*nrow(ldf),6)
   p = ggplot(data=cdf, aes(x=Sample, y=value, fill=variable)) +
      geom_bar(stat="identity", colour="black", size = 0.5) + ylab("Proportion") + xlab("") + scale_fill_manual(name = "Marker",values=colors) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(labels=rownames(ldf)) +
      geom_text(aes(label=TotalCells,y=1),vjust=-0.2)
   print(p)
   dev.off()
   return(p)
}

# Function to plot positive vs negative proportions for each marker and cell type combination
plot_NonCellTypeMarker = function(cell_type_markers, non_cell_type_markers, ldf, OutFileRoot) {
   # Get colors for each cell type marker
   colors = antibody_colors[cell_type_markers,"color"]
   
   # For each cell type marker
   for (cell_marker in cell_type_markers) {
      # Skip if the cell marker is "otherCell"
      if (cell_marker == "otherCell") next
      
      # For each non-cell type marker (e.g., Ki67, CD68)
      for (other_marker in non_cell_type_markers) {
         
         # Find all columns that match this cell type and marker
         matching_cols = grep(paste0("^", cell_marker, ".*", other_marker, "[+-]"), colnames(ldf), value=TRUE)
         
         if (length(matching_cols) == 0) {
            print(paste("No data found for", cell_marker, "and", other_marker))
            next
         }
         
         # Sum positive and negative counts
         positive_cols = grep(paste0(other_marker, "\\+"), matching_cols, value=TRUE)
         negative_cols = grep(paste0(other_marker, "\\-"), matching_cols, value=TRUE)
         
         # Create dataframe for plotting
         cdf = data.frame(
            Negative = rowSums(ldf[, negative_cols, drop=FALSE]),
            Positive = rowSums(ldf[, positive_cols, drop=FALSE])
         )
         
         # Calculate total cells for this cell type
         total = cdf$Negative + cdf$Positive
         
         # Skip if no cells of this type
         if (all(total == 0)) {
            print(paste("No cells found for", cell_marker))
            next
         }
         
         # Convert to proportions
         cdf = cdf/total
         cdf$Sample = rownames(cdf)
         
         # Write proportions to output file
         write.table(cdf, 
                    file = paste0(MainDir, OutFileRoot, gsub("/","-",cell_marker), "_", gsub("/","-",other_marker), ".txt"), 
                    quote = F, col.names = T, row.names = T, sep = "\t")
         
         # Reshape data for plotting
         cdf = melt(cdf, id="Sample")
         cdf = merge(cdf, 
                    data.frame(Sample = rownames(ldf), 
                             TotalCells = total,
                             PatientLymphomoidName = ldf$PatientLymphomoidName,
                             row.names = NULL), 
                    by = "Sample")
         
         cdf$variable = factor(cdf$variable, levels = c("Negative", "Positive"))
         cdf$Sample = factor(cdf$Sample, levels = rownames(ldf))
         cdf$TotalCells = paste0("N=", round(cdf$TotalCells))
         
         # Create stacked bar plot
         pdf(paste0(MainDir, OutFileRoot, gsub("/","-",cell_marker), "_", gsub("/","-",other_marker), ".pdf"), 
             2*nrow(ldf), 6)
         
         p = ggplot(data=cdf, aes(x=Sample, y=value, fill=variable)) +
            geom_bar(stat="identity", colour="black", size = 0.5) + 
            ylab("Proportion") + 
            xlab("") + 
            ggtitle(paste(cell_marker, other_marker, "status")) +
            scale_fill_manual(name = paste(other_marker, "status"),
                            values=c(colors[which(cell_marker==cell_type_markers)], "deepskyblue")) + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
            scale_x_discrete(labels=cdf$PatientLymphomoidName) +  # Use PatientLymphomoidName for labels
            geom_text(aes(label=TotalCells,y=1), vjust=-0.2)
         
         print(p)
         dev.off()
      }
   }
   return(p)
}

###############################

########## Main ##########
channels = read.table( file = ConfigTable, sep = "\t", header = T, quote = '' ,stringsAsFactors = F)
dir.create(paste0(MainDir,"Digital_IF_images/"), showWarnings = F)
dir.create(paste0(MainDir,"Classified_cells_tables/"), showWarnings = F)
if (Lymphomoids_to_process=="all"){ Lymphomoids_to_process = list.files(paste0( MainDir,"Lymphomoid_boundaries/" ), pattern = "Boundary.txt") }

# Set antibody colors based on whether mouse or human antibodies are used
if (grepl("mouse_channels.txt", ConfigTable)) {
  rownames(antibody_colors) = antibody_colors$antibody_mouse
} else if (grepl("human_channels.txt", ConfigTable)) {
  rownames(antibody_colors) = antibody_colors$antibody_human
}

# Extract cell type markers
cell_type_markerz = c("otherCell", channels$Antibody[channels$Is_CellType_Marker == TRUE])
other_markerz = channels$Antibody[channels$Is_CellType_Marker == FALSE]

# Create column names for summary table
# ToDo: add non-cell type markers in every possible negative and positive combination
all_colz = c("ImageName","PatientLymphomoidName", cell_type_markerz, "TotalCells")

# Initialize empty summary dataframe
tdf = data.frame(matrix(nrow = 0, ncol = length(all_colz), dimnames = list(NULL,all_colz)))

# Fix F4/80 naming inconsistency
colnames(tdf)[colnames(tdf)=="F4.80"] = "F4/80"

# Initialize vectors to track processed files for logging
quant_file_vec = c() # for log file
calib_file_vec = c() # for log file

# Process each lymphoid region
for (ll in Lymphomoids_to_process)
{
   cat("\n","Processing",ll,"...","\n" )
   
   ## Extract image and lymphoid region names from filename
   ImageName = sub("_[^_]+$", "", substr(ll,1,nchar(ll)-13))
   PatientLymphomoidName = substr(ll,nchar(ImageName)+2,nchar(ll)-13)
   
   # Load quantification data
   quant_file = paste0(MainDir,"Quantification/",ImageName,"/quantification/mesmer-",ImageName,"_merged.csv")
   quant_file_vec = c(quant_file_vec,substr(quant_file,nchar(MainDir)+1,nchar(quant_file)))
   if (!file.exists(quant_file)) { 
      cat("\n","Attention!",quant_file,"does not exist. Moving on to the next lymphomoid...","\n" )
      next }
   quant = read.csv(quant_file,stringsAsFactors = F)
   
   # Remove unnecessary morphological measurements from quantification data
   quant = quant[,!(colnames(quant) %in% c( "Area_nucleus","MajorAxisLength_nucleus","MinorAxisLength_nucleus","Eccentricity_nucleus","Solidity_nucleus","Extent_nucleus","Orientation_nucleus","X_centroid_cytoplasm","Y_centroid_cytoplasm","Area_cytoplasm","MajorAxisLength_cytoplasm","MinorAxisLength_cytoplasm","Eccentricity_cytoplasm","Solidity_cytoplasm","Extent_cytoplasm","Orientation_cytoplasm" ))]
   colnames(quant)[colnames(quant)=="F4.80"] = "F4/80"
   
   # Load calibration thresholds
   calib_file = paste0(MainDir,"Calibrated_thresholds/",ImageName,"_AllThresholds.txt")
   calib_file_vec = c(calib_file_vec,substr(calib_file,nchar(MainDir)+1,nchar(calib_file)))
   if (!file.exists(calib_file)) { 
      calib_file = paste0(MainDir,"Calibrated_thresholds/",ImageName,"_",PatientLymphomoidName,"_AllThresholds.txt")
      if (!file.exists(calib_file)) { 
      cat("\n","Attention! Calibration file not found. Moving on to the next lymphomoid...","\n" )
      next }}
   calib = read.table( file = calib_file, sep = " ", header = F, quote = '' ,stringsAsFactors = F)
   rownames(calib) = calib$V1
   
   # Load and parse lymphoid boundary coordinates
   pl = read.table(file = paste0(MainDir,"Lymphomoid_boundaries/",ll), sep = "-")
   pixelSize = as.numeric(substr(as.character(pl[2,"V1"]),10,nchar(as.character(pl[2,"V1"]))))
   plp_df = parse_boundary(pl,pixelSize)

   ## Configure thresholds for this lymphoid region
   ll_config = channels
   rownames(ll_config) = paste0(sub(" .*", "", ll_config$Channel_name),"_thresh")
   ll_config$calibrated_thresholds = calib[rownames(ll_config),"V2"]
   rownames(ll_config) = ll_config$Antibody

   ## Classify cells and identify those within lymphoid boundary
   quant = classify_cells(quant, ll_config)
   inPolygon = point.in.polygon(point.x = quant$spatial_1, point.y = quant$spatial_2, pol.x = plp_df$x, pol.y = plp_df$y)
   quant$in_lymphomoid = inPolygon==1

   ## Generate visualization plots if requested
   if (plot_IF_images)
   {
      fileName = paste0(MainDir,"Digital_IF_images/","IFimage_",ImageName,"_",PatientLymphomoidName,"_all.pdf")
      p = plot_digital_image(fileName, quant, plp_df, withOtherCells = T, onlyInLymphomoid = F, nonCellTypeMarker = NULL)
      fileName = paste0(MainDir,"Digital_IF_images/","IFimage_",ImageName,"_",PatientLymphomoidName,"_NoOtherCells_OnlyInLymphomoid.pdf")
      p = plot_digital_image(fileName, quant, plp_df, withOtherCells = F, onlyInLymphomoid = T, nonCellTypeMarker = NULL)

      for (marker in other_markerz)
      {
         if (marker != "DAPI") {
            fileName = paste0(MainDir,"Digital_IF_images/","IFimage_",ImageName,"_",PatientLymphomoidName,"_NoOtherCells_OnlyInLymphomoid_",gsub("/","-",marker),".pdf")
            p = plot_digital_image(fileName, quant, plp_df, withOtherCells = F, onlyInLymphomoid = T, nonCellTypeMarker = marker)
         }
      }
   }
   
   ## Clean up and save cell data
   quant$X_centroid_nucleus = NULL
   quant$Y_centroid_nucleus = NULL
   colnames(quant)[colnames(quant) %in% c( "spatial_1","spatial_2" )] = c( "x_centroid_nucleus_um","y_centroid_nucleus_um" )
   quant = quant[quant$in_lymphomoid,]
   
   # Skip if too few cells found
   if (nrow(quant) < 10 ) { 
      cat( "\n","*** Warning: less than 10 cells inside the boundary! check the IF images under 'Digital_IF_images/' to see what is the problem", "\n" )
      cat( "\n","*** Skipping",ll, "\n" )
      next
   }
   write.table(quant, file = paste0(MainDir,"Classified_cells_tables/Table_",ImageName,"_",PatientLymphomoidName,"_AllCells.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

   ## Create summary statistics for this lymphoid region
   ttdf = data.frame(matrix(0,nrow = 1, ncol = length(all_colz), dimnames = list(paste0(ImageName,"_",PatientLymphomoidName),all_colz)))
   colnames(ttdf)[colnames(ttdf)=="F4.80"] = "F4/80"
   ttdf[, c( "ImageName","PatientLymphomoidName" )] = c( ImageName,PatientLymphomoidName )

   # Add cell type counts
   ttdf[,names(table(quant$CellType_antibody))] = as.numeric(table(quant$CellType_antibody))
   ttdf[,"TotalCells"] = sum(as.numeric(table(quant$CellType_antibody)))

   # Get markers that aren't used for cell type classification
   non_celltype_markers = rownames(ll_config)[!ll_config$Is_CellType_Marker & ll_config$Antibody != "DAPI"]
   marker_names = ll_config[non_celltype_markers, "Antibody"]

   # Extract the counts for all combinations of cell_type marker and non-cell type markers
   for (cell_type_marker in cell_type_markerz) {
     
     # Generate all combinations of non-cell type markers
     combinations = expand.grid(lapply(marker_names, function(x) c(TRUE, FALSE)))
     colnames(combinations) = marker_names
     
     # Add columns for each combination to the ttdf data frame
     for (i in 1:nrow(combinations)) {
       combination = combinations[i, ]
       names(combination) = colnames(combinations) # necessary when there is only one marker
       combination_name = paste0(cell_type_marker, "_", paste0(names(combination), ifelse(combination, "+", "-"), collapse = "_"))
       
       # Count cells matching the combination, treating NA as FALSE
       ttdf[, combination_name] = sum(quant$CellType_antibody == cell_type_marker & 
         apply(as.data.frame(quant[, paste0(marker_names, "+")]), 1, function(row) { # as.data.frame handles the case in which there is only one marker (which would be otherwise interpreted as vector)
           all(ifelse(is.na(row), FALSE, row) == combination)
         }))
     }
   }

   tdf = rbind(tdf, ttdf)
}

## Saving and plotting at the image X lymphomoid level
dir.create(paste0(MainDir,"Results_ImageXLymphomoid_level/"),showWarnings = F)
save(tdf, file = paste0(MainDir,"Results_ImageXLymphomoid_level/SummaryTable_ImageXLymphomoid_AllCells.RData"))
write.table(tdf, file = paste0(MainDir,"Results_ImageXLymphomoid_level/SummaryTable_ImageXLymphomoid_AllCells.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
p = plot_CellTypeProportions(markers = cell_type_markerz, ldf = tdf, OutFileRoot = "Results_ImageXLymphomoid_level/ImageXLymphomoidLevel_CellTypeProportions_StackedBarplot")
p = plot_CellTypeProportions(markers = cell_type_markerz[cell_type_markerz!="otherCell"], ldf = tdf, OutFileRoot = "Results_ImageXLymphomoid_level/ImageXLymphomoidLevel_CellTypeProportions_ExclOtherCells_StackedBarplot")
p = plot_NonCellTypeMarker(
   cell_type_markers = cell_type_markerz[cell_type_markerz!="otherCell"], 
   non_cell_type_markers = other_markerz[other_markerz != "DAPI"],
   ldf = tdf, 
   OutFileRoot = "Results_ImageXLymphomoid_level/ImageXLymphomoidLevel_MarkerStatus_"
)

# ## Consolidating the same lymphomoid across images, saving and plotting
tdf$ImageName = NULL
ldf = aggregate(.~PatientLymphomoidName, data = tdf, FUN=mean)
rownames(ldf) = ldf$PatientLymphomoidName
save(ldf, file = paste0(MainDir,"SummaryTable_LymphomoidLevel_AllCells.RData"))
write.table(ldf, file = paste0(MainDir,"SummaryTable_LymphomoidLevel_AllCells.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
p = plot_CellTypeProportions(markers = cell_type_markerz, ldf = ldf, OutFileRoot = "LymphomoidLevel_CellTypeProportions_StackedBarplot")
p = plot_CellTypeProportions(markers = cell_type_markerz[cell_type_markerz!="otherCell"], ldf = ldf, OutFileRoot = "LymphomoidLevel_CellTypeProportions_ExclOtherCells_StackedBarplot")
p = plot_NonCellTypeMarker(
   cell_type_markers = cell_type_markerz[cell_type_markerz!="otherCell"], 
   non_cell_type_markers = other_markerz[other_markerz != "DAPI"],
   ldf = ldf, 
   OutFileRoot = "LymphomoidLevel_MarkerStatus_"
)
# ## Creating and saving log table
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


