# EVOS_M7000_ImageAnalysis
Custom codes to process and analyse M7000 images with DAPI, Actin, and up to 2 Proteins of Interest (PoI). 

Instructions: 

1. Initial File Preparation and Organization 

    File Structure Setup: Execute the Python script 1.py to establish the required directory structure for the pipeline.

    Image Formatting: Run the 2.ijm macro file (ImageJ/Fiji) on the raw image files located within both the DAPI and Actin folders.

2. Image Segmentation using Cellpose 

These steps outline the manual segmentation and correction process within the Cellpose application for both Actin and DAPI channels.

    Actin Channel Processing:

        Launch the Cellpose software.

        Open the image files for the Actin channel.

        Apply the pre-trained segmentation model, Actin_224.

        Manually review and correct the generated masks as needed.

        Save the corrected masks (Ctrl+R or equivalent command).

        Proceed to the next image (D or equivalent command).

        Repeat steps 4 through 6 for all subsequent images in the Actin directory.

    DAPI Channel Processing:

        Open the image files for the DAPI channel.

        Apply the pre-trained segmentation model, DAPI_234.

        Manually review and correct the generated masks as needed.

        Save the corrected masks (Ctrl+R).

        Proceed to the next image (D).

        Repeat steps 3 through 5 for all subsequent images in the DAPI directory.

3. Configuration and Data Measurement ⚙️

A. Data Consolidation and Script Transfer

    Consolidate Segmentation Data: Move all generated .zip files (containing the segmentation masks) from the Actin and DAPI output folders into their respective PoI (Points of Interest) analysis folders.

    Transfer Scripts and Configuration: Copy and paste the following required files into each PoI folder:

        conditions.txt

        4.py, 5.py, and 6.py

B. Configure conditions.txt

Open and modify the conditions.txt file to accurately map your experimental conditions to the image file prefixes. Adhere strictly to the following formatting rules:

    Structure: Each condition must be defined on a new line.

    Key-Value Delimitation: Use a colon (:) to separate the condition name (key) from its associated starter identifiers (values).

    Field Termination: Use a semicolon (;) to terminate the line.

    Starter Delimitation: Use commas (,) to separate multiple starter identifiers.

    Crucial Formatting Rule: DO NOT use any spaces immediately following the colon (:) delimiter.

    Example: 
    Arbitrary condition name 1:starter1,starter2;
    Arbitrary condition name 2:starter3,starter4,starter5;

C. Image Parameter Measurement

    Measure Parameters: Run the ImageJ/Fiji macro 3.ijm on the consolidated data to extract and measure the required morphological parameters from the segmented images.

4. Final Data Processing and Analysis 

    Process and Plot Data: Execute the final analysis scripts sequentially in the following order: 4.py, 5.py, and 6.py.
