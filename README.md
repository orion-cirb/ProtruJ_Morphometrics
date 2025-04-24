# ProtruJ_Morphometrics

* **Developed by:** Thomas & Héloïse
* **Developed for:** Laurent & Emmanuelle
* **Team:** Germain
* **Date:** April 2025
* **Software:** Fiji


### Images description

3 channels images, Nucleus, Protrusions, Bright field

A specific template image (patchTemplate.tif) available on GitHub 

     
### Macro description

1. Find patches based on the bright field channel.
2. Segment the Protrusions channel, save the result as "fluoImage.tif"  
3. Enlarge each patches Roi, keep only patch+protrusions, filter out small protrusions. 
4. Skeletonize the protrusions, run Geodesic map and run Local Thickness map and save the image for each pore.
5. Find the closest corner and also the angle based on the Centroid of the protrusion.
6. Write a global result file. 

### Dependencies

*patchTemplate.tif* files

*Multi-Template-Matching* + *IJ-OpenCV* FIJI Plug-in
*IJPB* FIJI Plug-in
### Version history

Version 0.9 released on April 24, 2025.
