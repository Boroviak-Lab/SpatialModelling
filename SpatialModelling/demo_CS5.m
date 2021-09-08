addpath(genpath('./Functions'))

%To do: function to transform input data same as for scaffold
%Precalculate cLim in GPInfer model. Return as .cLim 
%
%Load in the 3D scaffold
[OBJ1,section] = LoadCS5('3D');

%Plot the 3D model.
%If there are 2 inputs it will plot the surfaces of a specified
%set of tissues. Can specify figure number with third parameter if you have
%plots already open
tissues = 'all';
PlotEmbryoCS5(OBJ1,tissues);

%tissues = {'EmDisc','Am','ExMes','Tb','SYS','VE'};
tissues = {'all'} %,'Tb','ExMes','SYS','VE'}
PlotEmbryoCS5(OBJ1,tissues);
%PlotEmbryoCS5(OBJ1,{'Am','EmDisc'});
%PlotEmbryoCS5(OBJ1,tissues,1); %Specify figure no.
%view([0 90])
delete(findall(gcf,'Type','light'))
%camlight('left')
material shiny

%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
%Process the shots onto scaffold
[Output] = loadCS5Scaffold(D,Locations,Shots);

%Using PlotEmbryo funciton with a 4th/5th input projects shots on
%to the embryo (with some transparency)
PlotEmbryoCS5(OBJ1,{'all'},1,Output,'Location');
PlotEmbryoCS5(OBJ1,{'EmDisc'},2,Output,'Location');

%Now we optimise a GP model for all tissues (this generates the
%hyperparameters). This is the core of the inference, we must edit this to
%change the model e.g., different covariance functions. Different joint
%sharing of GPs over tissues etc.
[Output] = MarmosetGP_CS5(D,Output,'SOX2');
%...or training expression data
PlotEmbryoCS5(OBJ1,{'all'},1,Output,'Expression');
%Next we infer the GP mean function and project on to the surface of the 3D
%model
[Output] = MarmosetGPInfer_CS5(Output,OBJ1);
%Note that if we run MarmosetGP_CS5 we must always follow up with
%MarmosetGPInfer_CS5, otherwise the interpolated values will either be
%missing or from the gene we previiously run

%Now we can plot
PlotEmbryoCS5GP(Output,OBJ1,{'all'},2);
camlight('left')
PlotEmbryoCS5GP(Output,OBJ1,{'EmDisc'},3);
camlight('left')

%We have a helper function to load the default view for all tissues
[OBJ1b,a1,b1] = transformCS5(OBJ1,'all');
h = PlotEmbryoCS5GP(Output,OBJ1b,{'all'},2);
view(a1,b1)
camlight('left')

%print(['Plots/Demo_SOX2_CS5.pdf'],'-dpdf','-r1000');
%dpng and quality parameters

%Default view for AP axis
[OBJ2b,a2,b2] = transformCS5(OBJ1,'notall')
h = PlotEmbryoCS5GP(Output,OBJ2b,{'EmDisc','VE'},1);
view([a2,b2])
camlight('left')

%We can try now virtual IFs by loading a section. Here we load diagonal 2
%Options: 203,204,198,208,D1,D2
[OBJ1,section] = LoadCS5('D2');
%We have already optimised the hyperparameters for SOX2 so we simply infer
%values over the section
[Output] = MarmosetGPInfer_CS5(Output,OBJ1);
%Can use the plot function to see the section
PlotEmbryoCS5(OBJ1,'all');
%Now call the virtualIF functioin to look at virtual cross section
virtualIF_CS5(Output,OBJ1,{'all'},2)

%We can also do line plots
%[Output] = MarmosetGP_CS5(D,Output,'SOX2');
load('../Data/SpatialData/CS5_EmDisc.mat')
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
h = plotAP(Output, 1, [1,1,1], 'EmDisc');
title('SOX2')

%Virtual micropattern
load('../Data/SpatialData/CS5_EmDisc.mat')
virtualMicropatternCS5(D,Output,Line,{'SOX2','CDX1','MIXL1','T'},1);



