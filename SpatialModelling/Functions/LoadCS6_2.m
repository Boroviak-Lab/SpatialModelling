function [OBJ2B,sections] = LoadCS6_2(sections)

if strcmp(sections,'3D')==1
load('../Data/SpatialData/CS6_EI_twin_30042_HIGHRES.mat')
elseif strcmp(sections,'220')==1
load('../Data/SpatialData/Cross sections/CS6 - early (twin)/CS6_EI_twin_Revision_220.mat')
OBJ2B = OBJ2B_1; %Fine
elseif strcmp(sections,'212')==1
load('../Data/SpatialData/Cross sections/CS6 - early (twin)/CS6_EI_twin_Revision_212.mat')
OBJ2B = OBJ2B_1;
elseif strcmp(sections,'206')==1
load('../Data/SpatialData/Cross sections/CS6 - early (twin)/CS6_EI_twin_Revision_206.mat')
OBJ2B = OBJ2B_1;
elseif strcmp(sections,'D1')==1
%load('../Data/SpatialData/Cross sections/CS6 - early (twin)/CS6_EI_twin_Revision_040521_Diagonal1.mat')
%OBJ2B = OBJ2B_1;
%OBJ2B.objects([1,5,9,13,17,21,25]) = OBJ2B_1.objects([25,1,5,9,13,17,21]);
%OBJ2B.objects([4,8,12,16,20,24,28]) = OBJ2B_1.objects([28,5,8,12,16,20,24]);
elseif strcmp(sections,'D2')==1
load('../Data/SpatialData/Cross sections/CS6 - early (twin)/CS6_EI_twin_Revision_040521_Diagonal3.mat')
OBJ2B = OBJ2B_1;
elseif strcmp(sections,'D3')==1
load('../Data/SpatialData/Cross sections/CS6 - early (twin)/CS6_EI_twin_Revision_040521_Diagonal3.mat')
OBJ2B = OBJ2B_1;
else
load('../Data/SpatialData/CS6_EI_twin_30042_HIGHRES.mat')
end

%Am 1 EmDisc 5 VE 9 SYS 13 Stalk 17 Tb 21 ExMes 25
%D1
%EmD 1 VE 5 SYS 9 Stalk 13 Tb 17 Exmes 21 Am 25