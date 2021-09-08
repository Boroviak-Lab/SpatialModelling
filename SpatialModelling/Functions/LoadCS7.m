function [OBJ3,sections] = LoadCS7(sections)

if strcmp(sections,'3D')==1
load('../Data/SpatialData/200319_CS7_build_full_stalkcut_FORDYLAN.mat')
elseif strcmp(sections,'C2')==1
OBJ3_1 = load('../Data/SpatialData/Cross sections/CS7/CS7_cross_revisions.mat')
OBJ3 = OBJ3_1.OBJ3;
OBJ3.objects([1,5,9,13,17,27,31,4,8,12,16,20,30,34]) = OBJ3_1.OBJ3.objects([21,17,25,9,13,1,5,24,20,28,12,16,4,8]);
%OBJ3.objects([4,8,12,16,20,30,34]) = OBJ3_1.OBJ3.objects([24,20,28,12,16,4,8]);
else
load('../Data/SpatialData/200319_CS7_build_full_stalkcut_FORDYLAN.mat')
end