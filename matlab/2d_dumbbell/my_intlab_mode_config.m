%The path of the codes for switch between verified computing and approximate computing.
addpath('./mode_switch_interface/');  
addpath("./verified_eig_estimation/");

%The path of the library of verified eigenvalue estimation for matrix.
%addpath('./verified_eig_estimation/'); 

%The path of INTLAB toolbox and initialization.
addpath('~/app/Intlab_V12/');
%addpath('~/app/Intlab_V9/');
startintlab;

global INTERVAL_MODE;

%INTERVAL_MODE=1; for rigorous computing based on interval arithmetic.
%INTERVAL_MODE=0; for approximate computing with rounding error inside.
INTERVAL_MODE=1;

