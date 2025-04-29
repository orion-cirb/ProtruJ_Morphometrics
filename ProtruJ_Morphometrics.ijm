/*
 * Description: A tool to detect patches and compute metrics on protrusions.
 * Developed for: Laurent & Emmanuelle, Germain's team
 * Authors: Thomas Caille & Héloïse Monnet @ ORION-CIRB 
 * Date: April 2025
 * Repository: https://github.com/orion-cirb/ProtruJ_Morphometrics
 * Dependencies: patchTemplate.tif file in input directory (to download from GitHub repository)
 * 				 IJ-OpenCV-plugins, Multi-Template-Matching, and IJPB-plugins plugins (to install via Update Sites)
*/


/////////////// GLOBAL VARIABLES ///////////////
protrusionMinLengthMicron = 8 ; 
patchDilationPixel = 10 ;
fieldOfViewPixel = 100 ;
///////////////////////////////////////////////


// Clear ROI manager and close all opened images
roiManager("reset");
close("*");

// Ask for input directory
inputDir = getDirectory("Please select a directory containing images to analyze");

// Check patchTemplate.tif file is in input directory
if(!File.exists(inputDir + "patchTemplate.tif")){
	exit("Please download patchTemplate.tif file and drop it into your images directory.");
}

// Get all files in input directory
imgList = getFileList(inputDir);

// Create a global output directory
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
outDir = inputDir+"Results-"+(month+1)+"-"+dayOfMonth+"_"+hour+"-"+minute+File.separator();
if (!File.isDirectory(outDir)) {
	File.makeDirectory(outDir);
}

// Create a results csv file
resultsFile = File.open(outDir+ "results.csv");
print(resultsFile,"Image name, Patch ID, Protusion ID, Nearest corner, Max distance from patch (μm), Centroid direction (degree), Mean thickness (μm)\n");
File.close(resultsFile);

// Loop through all files with .nd2 extension in input directory
for (i = 0; i < imgList.length; i++) {
	if (endsWith(imgList[i], ".nd2")) {
		print("Analyzing image " + imgList[i] + "...");
		
    	// Create a specific output directory based on image name
		imgName = replace(imgList[i], ".nd2", "");
    	imgOutDir = outDir+ imgName + File.separator();
    	if (!File.isDirectory(imgOutDir)) {
			File.makeDirectory(imgOutDir);
		}
		
    	// Open patchTemplate.tif file
    	run("Bio-Formats Importer", "open=["+inputDir+"patchTemplate.tif] color_mode=Default view=Hyperstack stack_order=XYCZT");
    	getDimensions(patchWidth, patchHeight, channels, slices, frames);
    	// Open image
    	run("Bio-Formats Importer", "open=["+inputDir+imgList[i]+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
    	// Z-project
    	run("Z Project...", "projection=[Average Intensity]");
    	close(imgList[i]);
    	// Let the user rotate the image (all 3 channels)
    	Stack.setChannel(3);
    	run("Rotate... ");
 	 	
  		// Hide images during macro execution
    	setBatchMode(true);
    	
    	// Split channels
    	run("Split Channels");
    	close("C1-AVG_"+imgList[i]);
    	
    	// DETECT PATCHES ON BRIGHTFIELD CHANNEL
    	selectImage("C3-AVG_"+imgList[i]);
    	// Preprocess image
    	run("Median...", "radius=6");
    	//TODO: add run("Set Measurements...")
    	List.setMeasurements;
    	meanBright = List.getValue("Mean"); 
    	run("Subtract...", "value="+meanBright); 
    	run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 autoscale saturate");
    	rename("brightImage");
    	// Use Template matching algorithm to look for similar structures in the image than the one on patchTemplate.tif
    	run("Template Matching Image", "template=[patchTemplate.tif] image=[brightImage] rotate=[] matching_method=[Normalised 0-mean cross-correlation] number_of_objects=50 score_threshold=0.40 maximal_overlap=0 add_roi");
    	
    	// SEGMENT PROTRUSIONS CHANNEL
    	selectImage("C2-AVG_"+imgList[i]);
    	getPixelSize(unit, pixelWidth, pixelHeight);
    	// Preprocess image
    	run("Set Measurements...", "median redirect=None decimal=2");
    	List.setMeasurements;
    	
    	medianFluo = List.getValue("Median");
    	run("Median...", "radius=3"); 
    	run("Subtract...", "value="+medianFluo); 
    	
    	// Threshold image with automated method
    	setAutoThreshold("Triangle dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		// Save segmentation result
		saveAs("tiff", imgOutDir+"fluoBinaryMask");
		rename("fluoMask");
    	
    	// Initiate progress bar
 	 	run("Text Window...", "name=[Progress] width=25 height=1 monospaced");
 	 	
    	// Loop through all ROIs = all detected patches
    	for (j = 0; j < roiManager("count"); j++) {
 			// Update progress bar
   	 		print("[Progress]", "\\Update:"+(j+1)+"/"+roiManager("count")+" ("+((j+1)*100)/roiManager("count")+"%)\n");
   	 		
   	 		// Compute intensity standard deviation on brightfield channel, to filter out wrongly detected patches
   	 		selectImage("brightImage");
   	 		getDimensions(width, height, channels, slices, frames); 
			roiManager("select", j);
			Roi.getCoordinates(xpoints, ypoints);
			run("Set Measurements...", "standard redirect=None decimal=2");
			List.setMeasurements;
			
			// Patches are very bright with darker outline --> high intensity standard deviation, allowing us to filter out wrongly detected patches
			// Ignore wrongly detected patches + patches wrongly oriented + patches at the edges and delete associated ROI
			if ((List.getValue("StdDev") < 10000) || (ypoints[0] < (fieldOfViewPixel + patchHeight/2)) || (xpoints[0] < (fieldOfViewPixel + patchWidth/2))|| (xpoints[1] > width-(fieldOfViewPixel + patchWidth/2))|| (ypoints[1] > height-(fieldOfViewPixel + patchHeight/2))) { // TODO: use global variable=80 (see TODO below) instead of 150 + filter out patches on right and down borders as well using imgSize + remove (ypoints[2]-ypoints[1] > 110)?
				roiManager("delete");
				j=j-1;
			} else {				
				// Crop binary mask around enlarged patch
				selectImage("fluoMask");
				roiManager("select", j);
				run("Enlarge...", "enlarge="+fieldOfViewPixel); 
				run("Duplicate...", " ");
				getDimensions(widthEnlarge, heightEnlarge, channels, slices, frames);
				// Define coordinates of 5 corners of the patch
				cornerX = newArray(widthEnlarge/2+58,widthEnlarge/2+33,widthEnlarge/2-51,widthEnlarge/2-51,widthEnlarge/2+33);
				cornerY = newArray(heightEnlarge/2,heightEnlarge/2-36,heightEnlarge/2-36,heightEnlarge/2+35,heightEnlarge/2+35);
				
				// Create a pentagon using these coordinates
				makePolygon(cornerX[0],cornerY[0],cornerX[1],cornerY[1],cornerX[2],cornerY[2],cornerX[3],cornerY[3],cornerX[4],cornerY[4]);
				// Fill it in white
				setColor(255);
				run("Fill");
				run("Select None");
				// Filter out small objects to keep only the patch with attached protusions
				run("Set Measurements...", "mean redirect=None decimal=2"); // TOO: remove?
				run("Analyze Particles...", "size=3000-Infinity show=Masks exclude");
				run("Invert LUT");
				rename("protrusions_mask");
				
				// Skeletonize patch with attached protusions
				run("Duplicate...", " ");
				run("Skeletonize");
				rename("protrusions_skel");
				// Clear enlarged pentagon, to keep protusions skeleton only
				makePolygon(cornerX[0],cornerY[0],cornerX[1],cornerY[1],cornerX[2],cornerY[2],cornerX[3],cornerY[3],cornerX[4],cornerY[4]);
				run("Enlarge...", "enlarge="+patchDilationPixel+" pixel");
				run("Clear", "slice");
				run("Select None");
				
				// Create marker for the geodesic distance map: fill pentagon in white and clear everything outside
				run("Duplicate...", "title=marker");
				makePolygon(cornerX[0],cornerY[0],cornerX[1],cornerY[1],cornerX[2],cornerY[2],cornerX[3],cornerY[3],cornerX[4],cornerY[4]);
				run("Enlarge...", "enlarge="+patchDilationPixel+" pixel");
				run("Clear Outside");
				run("Fill");
				// Compute geodesic map to later retrieve distance between protusions extremity and the patch				
				run("Geodesic Distance Map", "marker=marker mask=protrusions_mask distances=[Verwer (12,17,27,38)] output=[32 bits] normalize");
				run("Multiply...", "value=" + pixelWidth);
				run("Properties...", "channels=1 slices=1 frames=1 pixel_width="+pixelWidth+" pixel_height="+pixelHeight+" voxel_depth=1");
				// Filter out small protrusions and analyze remaining one 	
				run("Set Measurements...", "min centroid redirect=protrusions_mask-geoddist decimal=2");
				
				selectImage("protrusions_skel");
				run("Analyze Particles...", "size="+protrusionMinLengthMicron+"-Infinity display clear overlay");
				nProtrusions = nResults;
				
				if (nProtrusions == 0) {
					params = imgName+","+(j+1)+",NaN,NaN,NaN,NaN,NaN";
					File.append(params, outDir + "results.csv");
				} else {

					// Compute the local thickness of the whole patch with extension and clear the polygon patch to extract only extensions thickness 
					selectImage("protrusions_mask"); 	
					run("Restore Selection"); 
					run("Clear", "slice"); 
					run("Select None");
					run("Local Thickness (masked, calibrated, silent)");
					run("Set Measurements...", "area mean min centroid redirect=protrusions_mask_LocThk decimal=2"); 
					selectImage("protrusions_skel"); 
					run("Analyze Particles...", "size="+protrusionMinLengthMicron+"-Infinity display overlay");
					
	
					// Loop through all the lines of the results table to retrieve information for each protrusion
					for (k = 0; k < nProtrusions; k++) { 
						 	
						// Get protrusion centroid and max intensity on geodesic map = distance between protusion extremity and the patch	
						maxDist = getResult("Max", k);
						centroidX = getResult("X", k)/pixelWidth;
						centroidY= getResult("Y", k)/pixelHeight;
							
						// Retrieve closest patch corner from protrusion centroid
						closestCornerDistance = 1000;
						closestCorner = 0;
						for (l = 0; l < cornerY.length; l++){
							tempDistance = sqrt(pow(centroidX-cornerX[l],2) + pow(centroidY-cornerY[l],2));
							if (tempDistance < closestCornerDistance)  {
								closestCornerDistance = tempDistance;
								closestCorner = l+1;	
							}
						}
						
						// Compute angle between horizontal line and the line from the center of the patch to the centroid of the protrusion	
						sizeX = getWidth();
						sizeY = getHeight();
						a = sizeX/2;
						b = sqrt(pow(sizeX - centroidX, 2) + pow(sizeY/2 - centroidY, 2));
						c = sqrt(pow(centroidX - sizeX/2, 2) + pow(centroidY - sizeY/2, 2));
						angle = acos((a*a+c*c-b*b) / (2*a*c))  *180 / 3.14 ;
						if (centroidY > sizeY/2) angle = 360-angle ;

						// Get the Mean value of the local thickness for each extension 
						meanLocThick = getResult("Mean", k+nProtrusions); 
						
						// Write results
						params = imgName+","+(j+1)+","+(k+1)+","+closestCorner+","+maxDist+","+angle+","+meanLocThick; 
						File.append(params, outDir + "results.csv");
					}

					// Save images
					selectImage("protrusions_mask");
					saveAs("tiff", imgOutDir+"patch"+(j+1)+"_mask");
					selectImage("protrusions_skel");
					saveAs("tiff", imgOutDir+"patch"+(j+1)+"_skel");
					selectImage("protrusions_mask-geoddist");
					saveAs("tiff", imgOutDir+"patch"+(j+1)+"_distMap");
					selectImage("protrusions_mask_LocThk");
					saveAs("tiff", imgOutDir+"patch"+(j+1)+"_locThck");
				}
				
				// Close all opened windows
				close("fluoMask-1");
				close("marker");
				close("*patch*");
				close("*protrusions*"); 
				close("Results");
	 		}
		}
		
	    // Save patches ROIs in a zip file
	    roiManager("save", imgOutDir + "patches.zip");
	    
	    // Close all windows
	   	print("[Progress]", "\\Close");
	    roiManager("reset");
	   	close("*");
	   	
	   	setBatchMode(false);
    }
}

print("Analysis done!");
