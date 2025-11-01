setBatchMode(true);

dir = getDirectory("Choose a folder containing your TIFF images");

cleanedDir = dir + "cleaned_tifs" + File.separator;
File.makeDirectory(cleanedDir);

list = getFileList(dir);

for (i=0; i<list.length; i++) {
if (endsWith(list[i], ".tif") || endsWith(list[i], ".TIF") || endsWith(list[i], ".tiff") || endsWith(list[i], ".TIFF")) {
path = dir + list[i];

    open(path);
    
    cleanedPath = cleanedDir + list[i];
    
    saveAs("Tiff", cleanedPath);
    
    close();
    
    print("Processed: " + list[i]);
}

}

setBatchMode(false);

print("Finished processing all TIFF files.");