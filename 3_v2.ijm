setBatchMode(true);

inputFolder=getDirectory("Choose a Directory containing your images and ROIs");

outputFolder=inputFolder+"output_ijm"+File.separator;
File.makeDirectory(outputFolder);
print("Output folder created: " + outputFolder);

allFiles=getFileList(inputFolder);

tiffFiles=newArray(0);
cellRoiFiles=newArray(0);
nucleusRoiFiles=newArray(0);

for (i=0; i<allFiles.length; i++) {
if (endsWith(allFiles[i], ".TIF")) {
baseName = substring(allFiles[i], 0, lengthOf(allFiles[i]) - 6);

    cellRoiName = baseName + "d2_rois.zip";
    nucleusRoiName = baseName + "d0_rois.zip";

    if (File.exists(inputFolder + cellRoiName) && File.exists(inputFolder + nucleusRoiName)) {
        tiffFiles = Array.concat(tiffFiles, newArray(allFiles[i]));
        cellRoiFiles = Array.concat(cellRoiFiles, newArray(cellRoiName));
        nucleusRoiFiles = Array.concat(nucleusRoiFiles, newArray(nucleusRoiName));
    } else {
        print("Warning: Incomplete file set found for basename " + baseName);
        print("    Looking for: " + allFiles[i] + ", " + cellRoiName + ", " + nucleusRoiName);
    }
}

}

if (tiffFiles.length == 0) {
print("Error: No complete sets of TIF, Cell ROI, and Nucleus ROI files found. Exiting.");
setBatchMode(false);
exit();
}

run("Set Measurements...", "area mean standard perimeter integrated redirect=None decimal=3");

function writeResult(filePath, isHeader, type, index_cell, index_nuc, stats) {
}

for (i=0; i<tiffFiles.length; i++) {
print("\n--- Starting analysis for file " + (i+1) + "/" + tiffFiles.length + " ---");

inputPathTiff = inputFolder + tiffFiles[i];
inputPathCellRoi = inputFolder + cellRoiFiles[i];
inputPathNucleusRoi = inputFolder + nucleusRoiFiles[i];

run("Bio-Formats", "open=["+inputPathTiff+"] autoscale color_mode=Default");
imageName=getTitle();

var unit, pw, ph, pd, pixel_inch;
    
getPixelSize(unit, pw, ph, pd);
print(unit,pw,ph,pd);
pixel_inch=1/pw;
    
if (unit=="inches") {
    run("Set Scale...", "distance="+pixel_inch+" known=25400 unit=um");
    print("NOTICE: Converted 'inches' to 'um' globally.");
}

outputFolderField = outputFolder + imageName + File.separator;
File.makeDirectory(outputFolderField);

roiManager("reset");
roiManager("Open", inputPathNucleusRoi);
nucleusCount = roiManager("count");

if (nucleusCount == 0) {
    print("No Nucleus ROIs found. Skipping.");
    run("Close");
    roiManager("reset");
    continue;
}

setForegroundColor(255, 255, 255);
newImage("Glaszczka", "8-bit black", 2048, 1536, 1);
for (k = 0; k < nucleusCount; k++) {
    selectImage("Glaszczka");
    roiManager("Select", k);
    run("Enlarge...", "enlarge=-1");
    run("Fill", "slice");
}

selectImage("Glaszczka");

maskPath = inputFolder + "Mask.png";
saveAs("PNG", maskPath);
print("Saved global watershed mask to: " + maskPath);

run("Close");


//new nucl mask file for correct cytoplasms generation

setForegroundColor(255, 255, 255);
newImage("Glaszczka2", "8-bit black", 2048, 1536, 1);
for (k = 0; k < nucleusCount; k++) {
    selectImage("Glaszczka2");
    roiManager("Select", k);
    run("Fill", "slice");
}

selectImage("Glaszczka2");

maskPath2 = inputFolder + "Mask_2.png";
saveAs("PNG", maskPath2);
print("Saved global watershed mask to: " + maskPath2);

run("Close");

//end


roiManager("reset");

roiManager("Open", inputPathCellRoi);
cellCount = roiManager("count");

if (cellCount == 0) {
    print("No Cell ROIs found. Skipping.");
    run("Close");
    roiManager("reset");
    continue;
}

if (isOpen("Results")) {selectWindow("Results"); run("Close");}

for (j = 0; j < cellCount; j++) {
    print("\n--- Processing Cell ROI index: " + j + "/" + (cellCount - 1) + " ---");

    newRoiStart = roiManager("count");

    selectWindow(imageName);
    roiManager("Select", j);
    run("Create Mask");
    rename("Mask");

    open(maskPath);

    imageCalculator("AND create", "Mask", "Mask.png");

    selectWindow("Mask");
    run("Close");
    selectWindow("Mask.png");
    run("Close");

    selectWindow("Result of Mask");
    run("Analyze Particles...", "size=0-Infinity add");

    newRoiCount = roiManager("count") - newRoiStart;
    print("Found " + newRoiCount + " nuclei inside Cell " + j + ".");

    if (newRoiCount > 0) {

        run("Clear Results");

        for (k = 0; k < newRoiCount; k++) {
            currentRoiIndex = newRoiStart + k;

            selectWindow(imageName);
            run("Select None");
            roiManager("Select", currentRoiIndex);
            
            //new: enlarge nucleus to avoid artifact
            run("Enlarge...", "enlarge=1");
            //end of new

            run("Measure");
        }
        saveAs("Results", outputFolderField + "nuclei_for_cell_" + j + ".csv");
        if (isOpen("Results")) {selectWindow("Results"); run("Close");}
        print("Saved nuclei measurements for " + newRoiCount + " nuclei for cell " + j);

        selectWindow("Result of Mask");
        run("Close");

        selectImage(imageName);
        roiManager("Select", j);
        run("Clear Results");
        run("Measure");
        saveAs("Results", outputFolderField + "cell_" + j + ".csv");
        if (isOpen("Results")) {selectWindow("Results"); run("Close");}
        print("Saved cell measurement for cell " + j);

        selectWindow(imageName);
        roiManager("Select", j);
        run("Create Mask");
        rename("Mask");

        open(maskPath2);

        imageCalculator("Substract create", "Mask", "Mask_2.png");

        selectWindow("Mask");
        run("Close");
        selectWindow("Mask_2.png");
        run("Close");

        selectWindow("Result of Mask");
        run("Create Selection");
        roiManager("Add");

        selectWindow(imageName);
        lastRoiIndex = roiManager("count") - 1;
        print(lastRoiIndex);
        roiManager("Select", lastRoiIndex);

        run("Clear Results");
        run("Measure");
        saveAs("Results", outputFolderField + "cytoplasm_for_cell_" + j + ".csv");
        if (isOpen("Results")) {selectWindow("Results"); run("Close");}
        print("Saved cytoplasm measurement for cell " + j);

        roiManager("Select", lastRoiIndex);
        roiManager("Delete");

        tempRoiIndices = newArray(newRoiCount);
        for (k = 0; k < newRoiCount; k++) {
            tempRoiIndices[k] = newRoiStart + k;
        }
        roiManager("Select", tempRoiIndices);
        roiManager("Delete");

        selectWindow("Result of Mask");
        run("Close");

        selectWindow(imageName);
        run("Select None");
        run("Duplicate...", "title=mask_background");
        setAutoThreshold("Triangle dark no-reset");
        setOption("BlackBackground", true);
        run("Convert to Mask");
        run("Invert");
        run("Create Selection");
        roiManager("Add");

        selectWindow(imageName);
        lastRoiIndex_background = roiManager("count") - 1;
        print(lastRoiIndex_background);
        roiManager("Select", lastRoiIndex_background);

        run("Clear Results");
        run("Measure");
        saveAs("Results", outputFolderField + "background_for_cell_" + j + ".csv");
        if (isOpen("Results")) {selectWindow("Results"); run("Close");}
        print("Saved background measurement for cell " + j);

        roiManager("Select", lastRoiIndex_background);
        roiManager("Delete");

        selectWindow("mask_background");
        run("Close");


    } else {
        selectWindow("Result of Mask");
        run("Close");

        print("Skipping cell " + j + " as no nuclei were found.");
        continue;
    }

    selectImage(imageName);
    roiManager("Deselect");
}

selectWindow(imageName);
run("Close");
roiManager("reset");

}

setBatchMode(false);
print("\nAll files processed successfully.");