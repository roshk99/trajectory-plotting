# trajectory-planning
***This project takes raw trajectories and creates a canal surface around them. This will allow for trajectory reproduction given an initial point.***

*Run `canal_visualization_gui` to play with the canal surface options with a GUI or run `main.m` to see the intermediate steps*

###Canal Visualization Options
* Dataset
	* Choose one of 4 datasets. Dataset 0 has artificial data and the remainder have experimental data
* Type of Canal
	* Choose either circles or ellipses as cross-sections
* Boundary Method (only applicable for circlular cross-sections)
	* Choose either fit to plane (`fit_boundaries.m`) or orthonormal vector method (`fit_boundaries2.m`)
* Cross-Section Method
	* Choose to align the cross-sections with either the TNB frames or the cross-sections themselves
* Trajectories
	* Choose a subset of the available trajectories to use to create the canal surface
* Plot Type
	* Plot with either individual circles or the surface (surface takes more graphics and is much slower)

###Animation and Reproduction (in progress)
* The animation and reproduction portions of this project are in progress. Check back later!

###Prerequisites
* MATLAB
* Robotics Toolbox (Peter Corke)
* FFMPEG (Media Editing Software)
* Bulk Rename Utility (or similar)

###Setup
* Follow instructions from the Robotics Toolbox for setup
    * Add toolbox to the path
    * Run startup.m in rvctools folder
* Download FFMPEG
    * Make sure it is added to your path

###Running
* Run main.m from the scripts folder with the desired options
    * Dataset number, type of canal, boundary calculation method, plot type

###Images to Video (Use Bulk Rename Utility or similar)
* Put all the image files from a given dataset in the same folder (need to rename so that the filenames do not overlap, but they stay in numeric order)
* Rename the files so that the first one is 0000.png and so on
* Open a terminal and navigate to the folder with all the image files
* Run the following command: `ffmpeg -r 10 -i %04d.png out.avi` replacing __out__ with your video file name
