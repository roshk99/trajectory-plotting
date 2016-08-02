function main()
% -----------------------------------------------------------------------
% A function that analyzes a specific dataset given certain parameters
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------
clc; clear; close all;

%addpath('scripts'); addpath('data');

set_to_run = 4; %0-4
fit_type = 'ellipses'; %circles or ellipses
plot_method = 'circles'; %circles or surface
how_many = 1:4; %1:4 for sets 0-2and 5, 1:5 for sets 3-4

%More parameters in runFunction.m
runFunction(set_to_run, fit_type, plot_method, how_many, false);
end