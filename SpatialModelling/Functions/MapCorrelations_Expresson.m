addpath(genpath('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/code/gpml-matlab-v3.6-2015-07-07'))

%Load CS7
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,ShotsCS5] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,ShotsCS5);

% %Load CS6
% [OBJ2,section] = LoadCS6('3D');
% [D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');
% [OutputCS6] = loadCS6Scaffold(D,Locations,ShotsCS6);
% 
% %Load CS7
% [OBJ3,section] = LoadCS7('3D');
% [D,Locations,XYZ,CellType,ShotsCS7] = LoadShots('CS7');
% [OutputCS7] = loadCS7Scaffold(D,Locations,ShotsCS7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now the cyno data
%Load in the IDs and correlations
LOCMarm = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset_newest/Loc_marm.csv');
IDMarm = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset_newest/Ano_marm.csv');
IDcyno   = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset_newest/Ano_cyno.csv');
CynoData = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset_newest/Data_cyno_Int.csv');
NNs = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset_newest/NearestNeighbour.csv');

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');

%Now we get expression
[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,'MIXL1');
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);
PlotEmbryoCS5GP(OutputCS5,OBJ1b,{'all'},1,[2,2,1]);
view(a1,b1)
camlight('left')

PlotEmbryoCS5GP(OutputCS5,OBJ1c,{'VE','EmDisc'},1,[2,2,3]);
view([138.5913,5.2765])
camlight('left')
view([c1,d1])

%Project NNs on to marmoset CS5
[OutputCS5] = NNCS5(OutputCS5,NNs,CynoData,IDcyno,'MIXL1');
[OutputCS5] = MarmosetGP_CS5Cyno(D,OutputCS5,'MIXL1');
[OutputCS5] = MarmosetGPInfer_CS5Cyno(OutputCS5,OBJ1);

PlotEmbryoCS5GP(OutputCS5,OBJ1b,{'all'},1,[2,2,2]);
view(a1,b1)
camlight('left')

PlotEmbryoCS5GP(OutputCS5,OBJ1c,{'VE','EmDisc'},1,[2,2,4]);
view([138.5913,5.2765])
camlight('left')
view([c1,d1])