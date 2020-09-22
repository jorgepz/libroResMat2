%% FMTrusS - Force's Method 2D Truss Solver  
% The FMTrusS is a simple GNU-Octave script that allows to solve 2D truss
% analysis problems using the Force's Method. The current version
% was generated for educational purposes. NO WARRANTY is provided and
% all errors and/or suggestions are welcome https://www.fing.edu.uy/~jorgepz/ .
addpath( [ pwd '/sources'] ); 

%% 1- Input data
% Input parameters for geometry and material definition of all the trusses.
% Applied loads and fixed degrees of freedom (supports) are also set here.
FMTrusSInput

%% 2- Pre-process
% Computation and storage of length, angle, E and A of all truss elements.
% Static indetermination computation. Equilibrium matrix assembly and 
% virtual states definition.
FMTrusSPreprocess

%% 3- Process
% Equilibrium equations solved, flexibility matrix assembly and
% virtual forces determined.
FMTrusSProcess

%% 4- Output
% Plot of truss structure with reactions, normal forces and external loads.
FMTrusSPlots
