function [D,Locations,XYZ,CellType,Shots] = LoadShots(stage)

%Load in the expressioin Data, global for all the things
D = importdata('./Data/SpatialData/X.csv');
Locations = importdata('./Data/SpatialData/IDs_and_locations.csv');
if strcmp(stage,'CS5')==1
    Shots = importdata('./Data/SpatialData/210528_CS5_Proj_shots_CP.csv');
elseif strcmp(stage,'CS6')==1
    Shots = importdata('./Data/SpatialData/210528_Proj_shots_CS6_EMBII_CP.csv')
elseif strcmp(stage,'CS62')==1
    Shots = importdata('./Data/SpatialData/210528_CS6_EI_newemb_Proj_shots_CP.csv')
elseif strcmp(stage,'CS7')==1
    Shots = importdata('../Data/SpatialData/210527_Proj_shots_CS7_CP.csv')
else
    disp('Wrong stage selected')
end
XYZ = Shots.data(:,1:3);
CellType = Shots.textdata(2:end,2);
