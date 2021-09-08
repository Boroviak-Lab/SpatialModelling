addpath(genpath('./Functions'))
%Load CS7
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,ShotsCS5] = LoadShots('CS5');
%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');
%Load CS7
[OBJ3,section] = LoadCS7('3D');
[D,Locations,XYZ,CellType,ShotsCS7] = LoadShots('CS7');

%Load in the IDs for the marmoset correlation matrices
LOCMarm = importdata('../Data/SpatialData/Amnioid/CrossCorr.csv')
IDMarm = importdata('../Data/SpatialData/Amnioid/Idents.csv')
%Correlations
C = LOCMarm.data; 

%Proecss CS6 data
[OutputCS6] = loadCS6ScaffoldCorrAmn(C,LOCMarm,IDMarm,ShotsCS6);
[OutputCS6] = MarmosetGP_CS6Corr(OutputCS6);


[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');


opac = 0.1;

m1 = []; m2 = []; m3 = []; m8 = []; m9 = [];

for i = 1 :size(C,2)
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,i,OBJ2);
m1 = [m1,OutputCS6.m1];
m2 = [m2,OutputCS6.m2];
m3 = [m3,OutputCS6.m3]; 
m8 = [m3,OutputCS6.m8]; 
m9 = [m9,OutputCS6.m9]; 
end

h1=figure(1)
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m2(:,1),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m3(:,1),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m8(:,1),'FaceColor','interp','LineStyle','none');
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);
axis off
ax = gca(h1)
set(ax,'clim',[.65,.9])
print(['./Plots/ESC1'],'-dpdf','-r1000');
close all

h2=figure(2)
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m2(:,1),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m3(:,1),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m8(:,1),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2c.objects(20).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m1(:,1),'FaceColor','interp','LineStyle','none');
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);
axis off
ax = gca(h2)
set(ax,'clim',[.65,.9])
print(['./Plots/ESC2'],'-dpdf','-r1000');
close all


close all
h1=figure(1)
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m2(:,2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m3(:,2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m8(:,2),'FaceColor','interp','LineStyle','none');
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);
axis off
ax = gca(h1)
set(ax,'clim',[.65,.8])
print(['./Plots/Am1'],'-dpdf','-r1000');
close all

h2=figure(2)
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m2(:,2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m3(:,2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m8(:,2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2c.objects(20).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',m1(:,2),'FaceColor','interp','LineStyle','none');
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);
axis off
ax = gca(h2)
set(ax,'clim',[.65,.8])
print(['./Plots/Am2'],'-dpdf','-r1000');
close all


