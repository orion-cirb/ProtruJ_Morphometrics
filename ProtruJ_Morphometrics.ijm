/*
 * Description: A tool to detect patches and compute metrics on protrusion.
 * Developed for: Laurent & Emmanuelle, Germain's team
 * Author: Thomas Caille & Héloïse Monnet @ ORION-CIRB 
 * Date: April 2025
 * Repository: https://github.com/orion-cirb/ProtruJ_Morphometrics
 * Dependencies: patchTemplate files available on GitHub repository
*/

/////////////// GLOBAL VARIABLES ///////////////

filamentSize = 8
patchEnlarge = 10


///////////////////////////////////////////////
// Ask for input directory
inputDir = getDirectory("Please select a directory containing images to analyze");
// Clear the ROI manager and close all opened images
roiManager("reset");
close("*");

if(!File.exists(inputDir + "patchTemplate.tif")){
	exit("Please download the patchTemplate.tif file and drop it into your image folder");
}
// Create an output result folder based on actual minute,hour ...
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
outDir = inputDir+"Results-"+(month+1)+"-"+dayOfMonth+"_"+hour+"-"+minute+File.separator();
if (!File.isDirectory(outDir)) {
	File.makeDirectory(outDir);
}

// Get all files in the input directory
list = getFileList(inputDir);

fileResults = File.open(outDir+ "results.csv");
		print(fileResults,"Image_Name, Pore_ID , Protusion_ID , Nearest_Corner , Distance_from_Pore (μm) , Protusion_Angle , Thickness_Mean (μm)\n");
		File.close(fileResults);



/////////////// Open and process images ///////////////
//////////////////////////////////////////////////////



// Loop through all files with .nd2 extension in the input directory
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], ".nd2")) {
    	
		imgName = replace(list[i],".nd2", "");
    	// Create an specific output directory based on the actual image name
    	imageOutDir = outDir+ imgName + File.separator();
    	if (!File.isDirectory(imageOutDir)) {
			File.makeDirectory(imageOutDir);
		}  	
    	// Open images, keep in mind, you must have the Template image in your input folder
    	 
    	run("Bio-Formats Importer", "open=["+inputDir + "patchTemplate.tif] color_mode=Default view=Hyperstack stack_order=XYCZT");
    	run("Bio-Formats Importer", "open=["+inputDir + list[i]+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
    	// Set the diplayed channel to 3 (bright)
    	run("Z Project...", "projection=[Average Intensity]");
    	Stack.setChannel(3);
    	// Let the user rotate the image (all 3 channels)
    	run("Rotate... ");
    	
    	// Initiate a progress bar
    	title = "[Progress]";
 	 	run("Text Window...", "name="+ title +" width=25 height=2 monospaced");
  		// Hide images during macro, greatly inprove speed
    	setBatchMode(true);
    	
    	// Found patches base on the 3rd channel (bright), to do so, smooth the image, apply a bandpass Filter to enhance the patches.
    	run("Split Channels");
    	selectImage("C3-AVG_"+list[i]);
    	run("Median...", "radius=6");
    	List.setMeasurements;
    	meanBright = List.getValue("Mean"); 
    	run("Subtract...", "value="+meanBright+""); 
    	run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 autoscale saturate");
    	rename("rawImage");
    	// Here we use the Template matching algorithm, we provide a template and it will look for similar structure in the image.
    	run("Template Matching Image", "template=[patchTemplate.tif] image=[rawImage] matching_method=[Normalised 0-mean cross-correlation] number_of_objects=50 score_threshold=0.40 maximal_overlap=0 add_roi show_result");
    	// Working on channel-2 (green)
    	selectImage("C2-AVG_"+list[i]);
    	getPixelSize(unit, pixelWidth, pixelHeight);
    	
    	// Preprocessing : enhance the Signal  substract the mean value of the image and apply a median filter
    	List.setMeasurements;
    	run("Median...", "radius=3");
    	medianFluo = List.getValue("Median");
    	run("Subtract...", "value="+medianFluo+""); 
    	
    	// Set an automatic threshold to segment the patches
    	setAutoThreshold("Triangle dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		
		//Save the resulting image 
		saveAs("tiff", imageOutDir + "fluoImage");
    	rename("fluoImage");
    	
    	
    	
/////////////// Create images of interest  ///////////////
/////////////////////////////////////////////////////////
    	
    	
    	
    	// Loop throught all the Roi's, all the patches detected
   	 	for (j = 0; j < roiManager("count"); j++) {
   	 		
 			// Progress bar
   	 		print(title, "\\Update:"+(j+1)+"/"+roiManager("count")+" ("+((j+1)*100)/roiManager("count")+"%)\n"+getBar((j+1), roiManager("count")));
   	 		// select the bright image and the first ROi, retrieve the Standard Deviation in the Roi's
   	 		selectImage("rawImage");
   	 		getDimensions(width, height, channels, slices, frames);
			roiManager("select",j);
			Roi.getCoordinates(xpoints, ypoints);
			List.setMeasurements;
			
			// As the patches are very bright and the outline darker, we filter patches by the difference between the max intensity and the min aka stdDev. 
			// We also filter out patches at the edges and delete associated ROi's 
			
			if ((List.getValue("StdDev") < 10000) || (ypoints[2]-ypoints[1] > 110) || (ypoints[0] < 150 ) || ( xpoints[0] < 150)) {
				roiManager("delete");
				j=j-1;
			} else {
				// define the coordinates of the 5 corners
				cornerX = newArray(245,220,136,136,220);
				cornerY = newArray(167,131,131,202,202);
				
				selectImage("fluoImage");
				run("Duplicate...", " ");
				roiManager("select",j);
				run("Enlarge...", "enlarge=80");
				run("Duplicate...", " ");
				
				// Create a pentagon covering the patch  
				makePolygon(cornerX[0],cornerY[0], cornerX[1],cornerY[1], cornerX[2],cornerY[2], cornerX[3],cornerY[3], cornerX[4],cornerY[4]);
				// Fill it in white
				setColor(255);
				run("Fill");
				run("Select None");
				// Filter object by size,in pixel, to keep only the mask of the patch with protusions
				run("Set Measurements...", "mean redirect=None decimal=2");
				run("Analyze Particles...", "size=+3000-Infinity show=Masks exclude ");
				run("Invert LUT");
				rename("ProtusionMask_of_Pore");
				run("Duplicate...", " ");
				
				// Performe the skeletonize of the patch with protusions
				run("Skeletonize");
				rename("Skel_of_Pore");
				makePolygon(cornerX[0],cornerY[0], cornerX[1],cornerY[1], cornerX[2],cornerY[2], cornerX[3],cornerY[3], cornerX[4],cornerY[4]);
				// Enlarge the "patch selection" and clear inside it, to keep only protusions skeleton,
				run("Enlarge...", "enlarge="+patchEnlarge+" pixel");
				run("Clear", "slice");
				run("Select None");
				
				// Create a pentagon covering the patch, clear all the signal outside the patch and then fill it in white, we will use it as a marker for the distance map 
				run("Duplicate...", " ");
				makePolygon(cornerX[0],cornerY[0], cornerX[1],cornerY[1], cornerX[2],cornerY[2], cornerX[3],cornerY[3], cornerX[4],cornerY[4]);
				run("Enlarge...", "enlarge="+patchEnlarge+" pixel");
				run("Clear Outside");
				run("Fill");
				rename("PoreMask_of_Pore");
				
				// Run the geodesic map to get the distance map of the protusions				
				run("Geodesic Distance Map", "marker=PoreMask_of_Pore mask=ProtusionMask_of_Pore distances=[Verwer (12,17,27,38)] output=[16 bits] normalize");
				// Analyse the particles inside the skeletonize image to retrieve Centroid and Max instensity (of the geodesic map)
				run("Set Measurements...", "mean redirect=None decimal=2");
				// Filter out the skeleton smaller than defined number of pixels
				selectImage("Skel_of_Pore");
				run("Analyze Particles...", "size="+filamentSize+"-Infinity display clear");
				
				// Store image size before the loop
				sizeX = getWidth();
				sizeY = getHeight();


/////////////// Compute and write results ///////////////
////////////////////////////////////////////////////////


				// Loop thought all the line of the result table to retrieve information on each extension
				for (k = 0; k < nResults; k++) {
					// Initialized the result variable
					params = ""+imgName;
					
					selectImage("Skel_of_Pore");
					run("Set Measurements...", "min centroid redirect=ProtusionMask_of_Pore-geoddist decimal=2");
					run("Properties...", "channels=1 slices=1 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");				
					run("Analyze Particles...", "size="+filamentSize+"-Infinity display clear overlay");
						
					// Define variables 
					Max = getResult("Max", k);
					X = getResult("X", k);
					Y = getResult("Y", k);
						
					// Loop and check all the corners and find the closest one
					closeCornerDistance = 500;
					for (l = 0; l <= (cornerY.length-1); l++){
						tempDistance = sqrt(pow((X - cornerX[l]),2) + pow((Y - cornerY[l]),2));
						if (tempDistance < closeCornerDistance)  {
							closeCornerDistance = tempDistance;
							closeCorner = l+1;	
						}
					}
					
					// Here we will compute the angle from the center of the patch to the centroid of the extension using the trigonometry circle 	
					a = sqrt(pow((sizeX/2 - sizeX),2)+ pow((sizeY/2 - sizeY/2),2));
					b = sqrt(pow((sizeX - X),2)+ pow((sizeY/2 - Y),2));
					c = sqrt(pow((X - sizeX/2),2) + pow((Y - sizeY/2),2));
					angle = (acos((((a*a)+(c*c))-(b*b))/(2*a*c)))*(180/3.14);
					if ( Y > sizeY/2) angle = 360-angle ;
						
					// Write results
					params += ","+j+1;
					params += ","+(k+1)+ ","+ closeCorner+","+(Max*pixelWidth)+ ","+angle;
					
					// Compute the local thickness of the whole patch with extension and clear the polygon patch to extract only extensions thickness 
					selectImage("ProtusionMask_of_Pore");
					run("Restore Selection");
					run("Clear", "slice");
					run("Select None");
					run("Local Thickness (masked, calibrated, silent)");
					run("Set Measurements...", "area mean min centroid median redirect=ProtusionMask_of_Pore_LocThk decimal=2");
					selectImage("Skel_of_Pore");
					run("Analyze Particles...", "size="+filamentSize+"-Infinity display clear overlay");
			
					// Get the Median value of the local thickness for each extension 
					Mean_locThick = getResult("Mean", k);							
					params += ","+ Mean_locThick;
					File.append(params, outDir + "results.csv");
					
					}
				// Save images 
				if (isOpen("ProtusionMask_of_Pore_LocThk")){
					selectImage("ProtusionMask_of_Pore_LocThk");
					saveAs("tiff", imageOutDir +(j+1)+ "_LocThk_Pore");
					selectImage("Skel_of_Pore");
					saveAs("tiff", imageOutDir +(j+1)+ "_Skel_Pore");
					selectImage("ProtusionMask_of_Pore-geoddist");
					saveAs("tiff", imageOutDir +(j+1)+ "_Distance_Pore");
				
				}
			// Close all the opened windows
			close("fluoImage-1");
			close("*of_Pore*");
			close("Results");
					   	
   	 		}
   	 	close("Progress");
    	}
    // Save the patches Roi's in a zip file
    roiManager("save", imageOutDir + "ROIs.zip");
    roiManager("reset");
   	print(title, "\\Close");
   	close("*");
   	setBatchMode(false);
    }
}   	



/////////////// Functions ///////////////
////////////////////////////////////////



function getBar(p1, p2) {
			        n = 20;
			        bar1 = "--------------------";
			        bar2 = "********************";
			        index = round(n*(p1/p2));
			        if (index<1) index = 1;
			        if (index>n-1) index = n-1;
			        return substring(bar2, 0, index) + substring(bar1, index+1, n);
  }
 