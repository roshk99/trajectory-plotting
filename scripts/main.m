function main()
% -----------------------------------------------------------------------
% A function that analyzes a specific dataset in the specified manner, see
% the parameters below
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------
close all;

set_to_run = 4; %0-3
fit_type = 'circles'; %circles or ellipses
plot_method = 'circles'; %circles or surface
how_many = 1:3;

run_func(set_to_run, fit_type, plot_method, how_many, false);

end