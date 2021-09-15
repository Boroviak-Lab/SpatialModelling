function [OBJ1,sections] = LoadCS5(sections)

if strcmp(sections,'3D')==1
load('./Data/SpatialData/CS5_revisions_HIGHRES.mat','OBJ1')
elseif strcmp(sections,'203')==1
load('./Data/SpatialData/Cross Sections/CS5/CS5_revisions_CROSS_203.mat')
OBJ1 = OBJ1_1;
elseif strcmp(sections,'204')==1
load('./Data/SpatialData/Cross Sections/CS5/CS5_revisions_CROSS_204.mat')
OBJ1 = OBJ1_1;
elseif strcmp(sections,'198')==1
load('./Data/SpatialData/Cross Sections/CS5/CS5_revisions_CROSS_198.mat')
OBJ1 = OBJ1_3;
elseif strcmp(sections,'208')==1
load('./Data/SpatialData/Cross Sections/CS5/CS5_revisions_CROSS_208.mat')
OBJ1 = OBJ1_4;
elseif strcmp(sections,'D1')==1
load('./Data/SpatialData/Cross Sections/CS5/CS5_revisions_DIAGONAL_1.mat')
OBJ1 = OBJ1_5;
elseif strcmp(sections,'D2')==1
load('./Data/SpatialData/Cross Sections/CS5/CS5_revisions_DIAGONAL_2.mat')
OBJ1 = OBJ1_6;
else
load('./Data/SpatialData/CS5_revisions_HIGHRES.mat','OBJ1')    
end
