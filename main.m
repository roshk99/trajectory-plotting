function main()
% -----------------------------------------------------------------------
% A function that analyzes a specific dataset in the specified manner, see
% the parameters below
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------
clc; clear; close all;

set_to_run = 4; %0-4
fit_type = 'circles'; %circles or ellipses
plot_method = 'surface'; %circles or surface
how_many = 1:5; %1:4 for sets 0-2, 1:5 for sets 3-4

runFunction(set_to_run, fit_type, plot_method, how_many, false);

end