*************************************************
Required libraries:
OpenCV
newmat
OpenMP (optional)

The code uses OpenMP for parallelization.   In order to turn this feature on, compile with the flag -fopenmp. If you do not wish to have parallelized code, simply compile without the flag and the #pragma omp instructions will be ignored. 

Libraries for linking:
gomp (optional)
newmat 
opencv_core
opencv_legacy
opencv_highgui
opencv_imgproc
opencv_calib3d

*************************************************************
************** Running and file formats *********************

The code loads grayscale images that represent the silhouette probability maps, and the camera calibration information.

Included in this folder is a folder named 'Sphere', with the following contents:
1. images in a format your version of OpenCV will open
2. cali.txt : This provides the paths of the images as well as the calibration information. My format is copied from the Middlebury Mview datasets (http://vision.middlebury.edu/mview/data/).

The format for each line is: "IMAGE_PATH\imagename.png k11 k12 k13 k21 k22 k23 k31 k32 k33 r11 r12 r13 r21 r22 r23 r31 r32 r33 t1 t2 t3 d0 d1 d2 d3 (d4) (d5) (d6) (d7)".    The projection matrix for that image is K*[R t]. The image origin is top-left, with x increasing horizontally, y vertically. Then, distortion coefficients, these can be either 4 or 8 elements; at least 4 are required.  If you already do the undistortion, just use zeros here (as I do in the Sphere folder).  I also use translation units as millimeters.
3. A folder named 'Experiment'.

In order to run the Sphere dataset with the code, you will have to change the path in cali.txt to the location on your machine.

Then, you will have to alter the Set.cpp file with the appropriate path of the 'Sphere' folder, so base_dir = path of your Sphere folder.

Then run without any arguments.  The results will be in the Experiment folder, and labeled according to the options specified in Set.cpp (more on that in a moment), and given the preliminary options I set this folder will be '6LS20_DS1_SP1'.  The reconstructions will be in the 'smoothed_files' folder, and in .ply format which can be viewed with the free Meshlab viewer.  The images in '6LS20_DS1_SP1' show the combination of the original images and the reconstruction images; white means that the pixel projects to a point of the reconstruction and that the original image registered the pixel as a silhouette pixel.  Magenta means that the original image registered the pixel as a silhouette pixel while the pixel does not project to a portion of the reconstruction. Blue means that the original image counts the pixel as not part of the silhouette, but that the pixel projects to a point on the reconstruction.

***************************************************************
**************** Changing options in Set.cpp ******************

Set.cpp currently looks like this:
switch (demo_number){
case 0: {
	base_dir =  "/home/jinglong/DemoData/Sphere";
	source_file = base_dir + "/cali.txt";

	write_directory = base_dir + "/Experiment/";

	downsample = 1.0;
	division = 20.0;

	pA = {-140, -120, -140};
	pB = {140, 120, 140};

		number_splits = 1;
	}	break;

} 

You can add your own directories and then call the program by: programname number_of_case

Now, an explanation of all of the options:

base_dir =  "/home/jinglong/DemoData/Sphere";
** Location of the folder containing the dataset**

source_file = base_dir + "/cali.txt";
** Name of the calibration file.  Note you can change this to whatever is convenient as long as the format is adhered to **

write_directory = base_dir + "/Experiment/";
** Where you are going to write the results. **

pA = {-140, -120, -140};
pB = {140, 120, 140};
** Specification of the bounding box; pA is the minimum x, y, z and pB the maximum x, y, z in mm. **

division = 20.0;
** This is the voxel resolution and in a cube-shaped voxel, would be the length of each side in millimeters.  We actually use spherical voxels, so this value is the diameter of the sphere. **

downsample = 1.0;
** It is possible to resize the image.  downsample = 1.0 means no resizing, downsample = 2.0 means that in each dimension the size if halved, etc.  This is sometimes helpful will really large images (5MB) and many images (90, for instance). **

number_splits = 1;
** This is the number of times of voxel size is reduced in the hierarchical approach in my PhD thesis.  Setting it to 1 means one division.  So given division = 20.0 in this example, the final voxel diameter will be 5mm.  If we set number_splits = 2, then the final voxel size is 2.5mm.  number_splits = 0 means that the voxel resolution will be 20mm. 


**************************************************************