# trajectory-planning
***This project takes raw position data of the end effector of the Kinova Jaco robot and analyzes it.***

** See main.m to reproduce the plots and videos in results folder **

###Data Analysis Steps
* Extracts raw position data, shifts, cuts, resizes, and smooths it
* Fits a canal surface over the data (a 3D envelope) using either circles or ellipses as cross-section and using one of two methods for calculating the surface
* For circles, creates a new trajectory based on initial points given within the first cross-section of the canal
* Animates a trajectory using the Jaco robot (inverse kinematics found using optimization so are not perfectly accurate)

###Plots generated 
* 4 data sets available, Set 0 has artificial data, and the Sets 1-3 have experimental data
* Raw and smoothed trajectories in 3D and components
* The canal surface surrounding the trajectories along with the mean trajectory (using either circles or ellipses and either method 1 or method 2 for boundary calculation)
* The generated trajectories on the canal surface plot
* An animation of the Jaco robot executing the trajectory

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
* Put all the image files from a given dataset in the same folder
(need to rename so that the filenames do not overlap, but they stay in 
numeric order)
* Rename the files so that the first one is 0000.png and so on
* Open a terminal and navigate to the folder with all the image files
* Run the following command: `ffmpeg -r 10 -i %04d.png out.avi` replacing __out__
with your video file name
