# trajectory-planning
***This project takes raw position data of the end effector of the Kinova Jaco robot and analyzes it.***

* Extracts position data, shifts, cuts, resizes, and smooths it
* Plots the components and 3D trajectories
* Simulates the movement of the robot along the trajectories (can convert to avi)
* Using optimization function for inverse kinematics of Jaco (it finds the elbow-down configuration and is not perfectly accurate)

**See results folder for plots and videos**

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
* Run main.m from the scripts folder
* 4 plots will appear and the script will create folders for the raw image files

###Images to Video (Use Bulk Rename Utility or similar)
* Put all the image files from a given dataset in the same folder
(need to rename so that the filenames do not overlap, but they stay in 
numeric order)
* Rename the files so that the first one is 0000.png and so on
* Open a terminal and navigate to the folder with all the image files
* Run the following command: `ffmpeg -r 10 -i %04d.png out.avi` replacing __out__
with your video file name
