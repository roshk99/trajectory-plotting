# trajectory-planning
***This project takes raw trajectories and creates a canal surface around them. This will allow for trajectory reproduction given an initial point.***

* Run `canal_visualization_gui` to play with the canal surface options with a GUI or run `main.m` to see the intermediate steps (parameters stored in `scripts/runFunction.m`)*

###Canal Visualization Options
* Dataset
	* Choose one of 5 datasets. Dataset 0 has artificial data and the remainder have experimental data
* Type of Canal
	* Choose either circles or ellipses as cross-sections
* Trajectories
	* Choose a subset of the available trajectories to use to create the canal surface
* Plot Type
	* Plot with either individual circles or the surface (surface takes more graphics and is much slower)

###Animation
* `scripts/runFunction.m` contains the demonstrations, mean trajectory, canal surface, and reproductions. This can be used along with the Robotics Toolbox for animation purposes. Obselete code to this end is in `scripts/animate.m`.

###Prerequisites
* MATLAB (along with Curve Fitting Toolbox)
* Robotics Toolbox (Peter Corke) for animation

###Running
* Run the following commands first:
	* `addpath('scripts'); addpath('data');`
* Run main.m from the scripts folder with the desired options
    * Dataset number, type of canal, boundary calculation method, plot type
* Adjust parameters in `scripts/runFunction.m`