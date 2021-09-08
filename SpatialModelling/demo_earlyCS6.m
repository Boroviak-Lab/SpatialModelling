addpath(genpath('./Functions'))

%Load in the 3D scaffold
[OBJ1,section] = LoadCS6_2('3D');

%Plot the 3D model.
%If there are 2 inputs it will plot the surfaces of a specified
%set of tissues. Can specify figure number with third parameter if you have
%plots already open
tissues = 'all';
PlotEmbryoCS6_2(OBJ1,tissues);
%PlotEmbryoCS5(OBJ1,tissues,1); %Specify figure no.
%camlight('right')

%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS62');
%Process the shots onto scaffold
[Output] = loadCS6Scaffold2(D,Locations,Shots);

[Output] = MarmosetGP_CS62(D,Output,'SOX2');
[Output] = MarmosetGPInfer_CS6(Output,OBJ1);

PlotEmbryoCS6GP_2(Output,OBJ1,{'Am','EmDisc','PGC'},1);

%Try some sections ...
%206,212,220,D2,D3 
[OBJ1,section] = LoadCS6_2('D3');
PlotEmbryoCS6_2(OBJ1,tissues,1);
[Output] = MarmosetGPInfer_CS6(Output,OBJ1);
virtualIF_CS62(Output,OBJ1,{'all'},1);



