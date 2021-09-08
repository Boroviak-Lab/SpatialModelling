%First plot gene expression ...
%Load CS5
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,Shots);
%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS6');
[OutputCS6] = loadCS6Scaffold(D,Locations,Shots);
%Load CS7
[OBJ3,section] = LoadCS7('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS7');
[OutputCS7] = loadCS7Scaffold(D,Locations,Shots);
OutputCS7.scalefactor = 600;
%...
[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');

OBJ1d = transformOBJ(OBJ1,pi/2,[0,1,0]);
OBJ1e = transformOBJ(OBJ1d,'flipX');
OBJ1f = transformOBJ(OBJ1e,'flipZ');
OBJ1test = transformOBJ(OBJ1f,6.0879,[0.6144,0.7715,0.1653]);
OBJ2c_part2 = transformCS6(OBJ2c,'flipY')
OBJ2test = transformOBJ(OBJ2c_part2,2.3663,[0.4229,0.7125,0.5599] );

[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,'T');
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);
[OutputCS6] = MarmosetGP_CS6_v3(D,OutputCS6,'T');
[OutputCS6] = MarmosetGPInfer_CS6_v3(OutputCS6,OBJ2);
[OutputCS7] = MarmosetGP_CS7_v3(D,OutputCS7,'T');
[OutputCS7] = MarmosetGPInfer_CS7_v3(OutputCS7,OBJ3);
m1 = OutputCS5.cLim;
m2 = OutputCS6.cLim;
m3 = OutputCS7.cLim;
ma1s = max(max([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));
mi1s = min(min([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));
sf = 1;
isf = 1/sf;
mi1 = min([m1,m2,m3]);
ma1 = max([m1,m2,m3]);
m2 = max([ma1,0.1*isf]);
m1 = max([mi1,0]);
mi2 = max([ma1s,0.1*isf]);
mi1 = max([mi1s,0]);
ma1s = max(max([OutputCS5.cLim2;OutputCS6.cLim2;OutputCS7.cLim2]));
mi1s = min(min([OutputCS5.cLim2;OutputCS6.cLim2;OutputCS7.cLim2]));


h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1test,{'EmDisc'},2,[5,8,17,18,25,26]);
view([375.7415,-26.7022])
camlight('left')
ax = gca(h2)
ax.CLim = [mi1 mi2]*sf;
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;
h4 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2c,{'EmDisc','Stalk'},2,[5,8,19,20,27,28]);
view([-70.8201,90])
camlight('left')
ax = gca(h4)
ax.CLim = [mi1 mi2]*sf;
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;
cb=colorbar;
cb.Position = cb.Position + [0.06 0.05 -0.002 -0.1];
cb.Color = 'k';
cb.FontSize = 12;
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(end)])
end

print(['./Plots/T'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/POU5F1.fig'])
close all


[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,'MIXL1');
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);
[OutputCS6] = MarmosetGP_CS6_v3(D,OutputCS6,'MIXL1');
[OutputCS6] = MarmosetGPInfer_CS6_v3(OutputCS6,OBJ2);
[OutputCS7] = MarmosetGP_CS7_v3(D,OutputCS7,'MIXL1');
[OutputCS7] = MarmosetGPInfer_CS7_v3(OutputCS7,OBJ3);
m1 = OutputCS5.cLim;
m2 = OutputCS6.cLim;
m3 = OutputCS7.cLim;
ma1s = max(max([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));
mi1s = min(min([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));
sf = 1;
isf = 1/sf;
mi1 = min([m1,m2,m3]);
ma1 = max([m1,m2,m3]);
m2 = max([ma1,0.1*isf]);
m1 = max([mi1,0]);
mi2 = max([ma1s,0.1*isf]);
mi1 = max([mi1s,0]);
ma1s = max(max([OutputCS5.cLim2;OutputCS6.cLim2;OutputCS7.cLim2]));
mi1s = min(min([OutputCS5.cLim2;OutputCS6.cLim2;OutputCS7.cLim2]));

h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1test,{'EmDisc'},2,[5,8,17,18,25,26]);
view([375.7415,-26.7022])
camlight('left')
ax = gca(h2)
ax.CLim = [mi1 mi2]*sf;
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;
h4 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2c,{'EmDisc','Stalk'},2,[5,8,19,20,27,28]);
view([-70.8201,90])
camlight('left')
ax = gca(h4)
ax.CLim = [mi1 mi2]*sf;
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;
cb=colorbar;
cb.Position = cb.Position + [0.06 0.05 -0.002 -0.1];
cb.Color = 'k';
cb.FontSize = 12;
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(end)])
end

print(['./Plots/MIXL1'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/POU5F1.fig'])
close all

[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,'SOX2');
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);
[OutputCS6] = MarmosetGP_CS6_v3(D,OutputCS6,'SOX2');
[OutputCS6] = MarmosetGPInfer_CS6_v3(OutputCS6,OBJ2);
[OutputCS7] = MarmosetGP_CS7_v3(D,OutputCS7,'SOX2');
[OutputCS7] = MarmosetGPInfer_CS7_v3(OutputCS7,OBJ3);
m1 = OutputCS5.cLim;
m2 = OutputCS6.cLim;
m3 = OutputCS7.cLim;
ma1s = max(max([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));
mi1s = min(min([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));
sf = 1;
isf = 1/sf;
mi1 = min([m1,m2,m3]);
ma1 = max([m1,m2,m3]);
m2 = max([ma1,0.1*isf]);
m1 = max([mi1,0]);
mi2 = max([ma1s,0.1*isf]);
mi1 = max([mi1s,0]);
ma1s = max(max([OutputCS5.cLim2;OutputCS6.cLim2;OutputCS7.cLim2]));
mi1s = min(min([OutputCS5.cLim2;OutputCS6.cLim2;OutputCS7.cLim2]));

h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1test,{'EmDisc'},2,[5,8,17,18,25,26]);
view([375.7415,-26.7022])
camlight('left')
ax = gca(h2)
ax.CLim = [mi1 mi2]*sf;
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;
h4 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2c,{'EmDisc','Stalk'},2,[5,8,19,20,27,28]);
view([-70.8201,90])
camlight('left')
ax = gca(h4)
ax.CLim = [mi1 mi2]*sf;
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;
cb=colorbar;
cb.Position = cb.Position + [0.06 0.05 -0.002 -0.1];
cb.Color = 'k';
cb.FontSize = 12;
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(end)])
end

print(['./Plots/SOX2'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/POU5F1.fig'])
close all

%Load in pre-existing runs to save time
load(['../Data/SpatialData/PrimedNiave/Pre_Marm_Int.mat'],'Pre')
load(['../Data/SpatialData/PrimedNiave/CS5Mapping_Marm_Int.mat'])
load(['../Data/SpatialData/PrimedNiave/CS6Mapping_Marm_Int.mat'])
IDs = importdata('../Data/SpatialData/PrimedNiave/InVit_Idents.csv')
Cls = importdata('../Data/SpatialData/PrimedNiave/InVit_Cl5.csv')

IDs = IDs.textdata(2:end,2);
indn = find(strcmp(IDs,'PLAXA')==1 | strcmp(IDs,'EMS3_PLAXA')==1);
cln = Cls.data(indn,1);

indp = find(strcmp(IDs,'Conventional')==1 | strcmp(IDs,'ESC_CMESC')==1 | strcmp(IDs,'ESC_conv2')==1);
clp = Cls.data(indp,1);
indp1 = indp(find(clp==3));
indp2 = indp(find(clp==8));

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
LOCMarm = importdata('../Data/SpatialData/PrimedNiave/Marm_LOCS.csv')
IDMarm = importdata('../Data/SpatialData/PrimedNiave/Marm_Idents.csv')
C = importdata('../Data/SpatialData/PrimedNiave/Marm_InVit_Int.csv');
%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');

OBJ1d = transformOBJ(OBJ1,pi/2,[0,1,0]);
OBJ1e = transformOBJ(OBJ1d,'flipX');
OBJ1f = transformOBJ(OBJ1e,'flipZ');

%Also load the blastocyst
OBJ0a=importdata('../Data/SpatialData/Blastocyst_halfcut_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ0b=importdata('../Data/SpatialData/Blastocyst_slice_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
meh = load('../Data/SpatialData/lBlHR.mat')
OBJ0b = meh.O4;
OBJ0b.vertices = Quaternion3(pi*1.28,[0,1,0],OBJ0b.vertices);
OBJ1test = transformOBJ(OBJ1f,6.0879,[0.6144,0.7715,0.1653]);
OBJ2c_part2 = transformCS6(OBJ2c,'flipY')
OBJ2test = transformOBJ(OBJ2c_part2,2.3663,[0.4229,0.7125,0.5599] );

load('../Data/SpatialData/eBlHR2.mat')

%Now begin plotting
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
%view([5.5713,-90])
camlight('left')
material shiny
colormap(parula);
%set(gca,'clim',[.65,.76] )

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indp),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indp),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indp),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indp),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indp),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indp),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indp),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indp),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp)))./(size(Pre.C5(:,indp),2)*size(Pre.C5(:,indp),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indp)))./(size(Pre.C3(:,indp),2)*size(Pre.C3(:,indp),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indp)))./(size(Pre.C1(:,indp),2)*size(Pre.C1(:,indp),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);

set(h1,'clim',[.45,.65] )
set(h2,'clim',[.45,.65] )
set(h3,'clim',[.45,.65] )
set(h4,'clim',[.45,.65] )
set(h5,'clim',[.45,.65] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Conv'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv.fig'])
close all

h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp)))./(size(Pre.C5(:,indp),2)*size(Pre.C5(:,indp),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indp)))./(size(Pre.C0(:,indp),2)*size(Pre.C0(:,indp),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h5,'clim',[.45,.65] )
print(['./Plots/Conv_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_icm.fig'])
close all 

h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny 
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view([5.5713,-90])
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indp1),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indp1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indp1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indp1),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);


h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indp1),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indp1),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indp1),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indp1),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp1)))./(size(Pre.C5(:,indp1),2)*size(Pre.C5(:,indp1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indp1)))./(size(Pre.C3(:,indp1),2)*size(Pre.C3(:,indp1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indp1)))./(size(Pre.C1(:,indp1),2)*size(Pre.C1(:,indp1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);

set(h1,'clim',[.45,.65] )
set(h2,'clim',[.45,.65] )
set(h3,'clim',[.45,.65] )
set(h4,'clim',[.45,.65] )
set(h5,'clim',[.45,.65] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Conv_1'],'-dpdf','-r1000');
savefig(['./Plots/Conv_1.fig'])
close all


h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp1)))./(size(Pre.C5(:,indp1),2)*size(Pre.C5(:,indp1),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indp1)))./(size(Pre.C0(:,indp1),2)*size(Pre.C0(:,indp1),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
%view([375.7415,-26.7022])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h5,'clim',[.45,.65] )
print(['./Plots/Conv_1_icm'],'-dpdf','-r1000');
savefig(['./Plots/Conv_1_icm.fig'])
close all

h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);


h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view([5.5713,-90])
view([-70.8201,90])

camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indp2),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indp2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indp2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indp2),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indp2),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indp2),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indp2),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indp2),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny 
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp2)))./(size(Pre.C5(:,indp2),2)*size(Pre.C5(:,indp2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indp2)))./(size(Pre.C3(:,indp2),2)*size(Pre.C3(:,indp2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indp2)))./(size(Pre.C1(:,indp2),2)*size(Pre.C1(:,indp2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);

set(h1,'clim',[.45,.65] )
set(h2,'clim',[.45,.65] )
set(h3,'clim',[.45,.65] )
set(h4,'clim',[.45,.65] )
set(h5,'clim',[.45,.65] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Conv_2'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_2.fig'])
close all


h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp2)))./(size(Pre.C5(:,indp2),2)*size(Pre.C5(:,indp2),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indp2)))./(size(Pre.C0(:,indp2),2)*size(Pre.C0(:,indp2),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h5,'clim',[.45,.65] )
print(['./Plots/Conv_2_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_2_icm.fig'])
close all

h = figure(3)
h1=subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny 
colormap(parula);

h2=subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices', OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view([5.5713,-90])
view([-70.8201,90])
camlight('left')
material shiny 
colormap(parula);

h3=subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indn),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indn),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indn),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indn),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4=subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indn),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indn),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indn),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indn),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny 
colormap(parula);

h5=subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn)))./(size(Pre.C5(:,indn),2)*size(Pre.C5(:,indn),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indn)))./(size(Pre.C3(:,indn),2)*size(Pre.C3(:,indn),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indn)))./(size(Pre.C1(:,indn),2)*size(Pre.C1(:,indn),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);

set(h1,'clim',[.45,.65] )
set(h2,'clim',[.45,.65] )
set(h3,'clim',[.45,.65] )
set(h4,'clim',[.45,.65] )
set(h5,'clim',[.45,.65] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Niave'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Niave.fig'])
close all


h = figure(3)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn)))./(size(Pre.C5(:,indn),2)*size(Pre.C5(:,indn),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indn)))./(size(Pre.C0(:,indn),2)*size(Pre.C0(:,indn),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h5,'clim',[.45,.65] )
print(['./Plots/Niave_cim'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Niave_icm.fig'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now we get on to human...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['../Data/SpatialData/PrimedNiave/Pre_PSC_Int.mat'],'Pre')
load(['../Data/SpatialData/PrimedNiave/CS5Mapping_HumanRun_Int.mat'])
load(['../Data/SpatialData/PrimedNiave/CS6Mapping_HumanRun_Int.mat'])
IDs = importdata('../Data/SpatialData/PrimedNiave/PSC_Idents.csv')
Cls = importdata('../Data/SpatialData/PrimedNiave/PSC_Cl5.csv')

IDs = IDs.textdata(2:end,2);
indn = find(strcmp(IDs,'Naive')==1 | strcmp(IDs,'EMS3_PLAXA')==1);
cln = Cls.data(indn,1);

indp = find(strcmp(IDs,'Primed')==1 | strcmp(IDs,'ESC_CMESC')==1 | strcmp(IDs,'ESC_conv2')==1);
clp = Cls.data(indp,1);

indp1 = indp(find(clp==0));
indp2 = indp(find(clp==3));
indp3 = indp(find(clp==12));

indn1 = indn(find(cln==6));
indn2 = indn(find(cln==4));
indn3 = indn(find(cln==2));
indn4 = indn(find(cln==12));

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
LOCMarm = importdata('../Data/SpatialData/PrimedNiave/Marm_LOCS.csv')
IDMarm = importdata('../Data/SpatialData/PrimedNiave/Marm_Idents.csv')
C = importdata('../Data/SpatialData/PrimedNiave/Marm_InVit_Int.csv');
%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');

OBJ1d = transformOBJ(OBJ1,pi/2,[0,1,0]);
OBJ1e = transformOBJ(OBJ1d,'flipX');
OBJ1f = transformOBJ(OBJ1e,'flipZ');

%Also load the blastocyst
OBJ0a=importdata('../Data/SpatialData/Blastocyst_halfcut_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ0b=importdata('../Data/SpatialData/Blastocyst_slice_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')

meh = load('../Data/SpatialData/lBlHR.mat')
OBJ0b = meh.O4;
OBJ0b.vertices = Quaternion3(pi*1.28,[0,1,0],OBJ0b.vertices);
OBJ1test = transformOBJ(OBJ1f,6.0879,[0.6144,0.7715,0.1653]);
OBJ2c_part2 = transformCS6(OBJ2c,'flipY')
OBJ2test = transformOBJ(OBJ2c_part2,2.3663,[0.4229,0.7125,0.5599] );

%Conventional 
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([5.5713,-90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indp),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indp),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indp),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indp),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indp),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indp),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indp),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indp),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indp),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny
colormap(parula);


%5=tb,3=hyp,1=epi
h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp)))./(size(Pre.C5(:,indp),2)*size(Pre.C5(:,indp),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indp)))./(size(Pre.C3(:,indp),2)*size(Pre.C3(:,indp),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indp)))./(size(Pre.C1(:,indp),2)*size(Pre.C1(:,indp),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);

set(h1,'clim',[.45,.68] )
set(h2,'clim',[.45,.68] )
set(h3,'clim',[.45,.68] )
set(h4,'clim',[.45,.68] )
set(h5,'clim',[.45,.68] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
filem = ['Conv.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/Conv_Human'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human.fig'])
close all

%Now the ICM
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp)))./(size(Pre.C5(:,indp),2)*size(Pre.C5(:,indp),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indp)))./(size(Pre.C0(:,indp),2)*size(Pre.C0(:,indp),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h5,'clim',[.45,.68] )
filem = ['Conv.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/Conv_Human_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human_icm.fig'])
close all

%Now look at subclusters%%%%%%%%%%%%%%%%%%%%%
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([5.5713,-90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indp1),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indp1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indp1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indp1),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indp1),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp1),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indp1),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indp1),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indp1),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indp1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny 
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp1)))./(size(Pre.C5(:,indp1),2)*size(Pre.C5(:,indp1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indp1)))./(size(Pre.C3(:,indp1),2)*size(Pre.C3(:,indp1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indp1)))./(size(Pre.C1(:,indp1),2)*size(Pre.C1(:,indp1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h1,'clim',[.45,.68] )
set(h2,'clim',[.45,.68] )
set(h3,'clim',[.45,.68] )
set(h4,'clim',[.45,.68] )
set(h5,'clim',[.45,.68] )
h1.XLim = h1.XLim*0.9;%
h1.YLim = h1.YLim*0.9;
print(['./Plots/Conv_Human_1'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human_1.fig'])
close all

%Now the ICM
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp1)))./(size(Pre.C5(:,indp1),2)*size(Pre.C5(:,indp1),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indp1)))./(size(Pre.C0(:,indp1),2)*size(Pre.C0(:,indp1),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h5,'clim',[.45,.68] )
filem = ['Conv.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/Conv_Human_1_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human_1_icm.fig'])
close all


%Next cluster
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indp2),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indp2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indp2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indp2),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indp2),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp2),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indp2),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indp2),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indp2),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indp2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny 
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp2)))./(size(Pre.C5(:,indp2),2)*size(Pre.C5(:,indp2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indp2)))./(size(Pre.C3(:,indp2),2)*size(Pre.C3(:,indp2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indp2)))./(size(Pre.C1(:,indp2),2)*size(Pre.C1(:,indp2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);

set(h1,'clim',[.45,.68] )
set(h2,'clim',[.45,.68] )
set(h3,'clim',[.45,.68] )
set(h4,'clim',[.45,.68] )
set(h5,'clim',[.45,.68] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Conv_Human_2'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human_2.fig'])
close all

%Now the ICM
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp2)))./(size(Pre.C5(:,indp2),2)*size(Pre.C5(:,indp2),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indp2)))./(size(Pre.C0(:,indp2),2)*size(Pre.C0(:,indp2),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h5,'clim',[.45,.68] )
print(['./Plots/Conv_Human_2_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human_2_icm.fig'])
close all

%%%%%%%%%%%%%%
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp3),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp3),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indp3),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indp3),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indp3),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indp3),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indp3),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indp3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indp3),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp3),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indp3),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indp3),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indp3),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indp3),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indp3),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indp3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny 
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp3)))./(size(Pre.C5(:,indp3),2)*size(Pre.C5(:,indp3),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indp3)))./(size(Pre.C3(:,indp3),2)*size(Pre.C3(:,indp3),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indp3)))./(size(Pre.C1(:,indp3),2)*size(Pre.C1(:,indp3),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);

set(h1,'clim',[.45,.68] )
set(h2,'clim',[.45,.68] )
set(h3,'clim',[.45,.68] )
set(h4,'clim',[.45,.68] )
set(h5,'clim',[.45,.68] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Conv_Human_3'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human_3.fig'])
close all

%Now the ICM
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indp3)))./(size(Pre.C5(:,indp3),2)*size(Pre.C5(:,indp3),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indp3)))./(size(Pre.C0(:,indp3),2)*size(Pre.C0(:,indp3),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h5,'clim',[.45,.68] )
filem = ['Conv.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/Conv_Human_3_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Conv_Human_3_icm.fig'])
close all

%Niave%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])

camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indn),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indn),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indn),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indn),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indn),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indn),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indn),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indn),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indn),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn)))./(size(Pre.C5(:,indn),2)*size(Pre.C5(:,indn),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indn)))./(size(Pre.C3(:,indn),2)*size(Pre.C3(:,indn),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indn)))./(size(Pre.C1(:,indn),2)*size(Pre.C1(:,indn),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);

set(h1,'clim',[.45,.75] )
set(h2,'clim',[.45,.75] )
set(h3,'clim',[.45,.75] )
set(h4,'clim',[.45,.75] )
set(h5,'clim',[.45,.75] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human'],'-dpdf','-r1000');
savefig(['./Plots/Naive_Human.fig'])


close all
h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn)))./(size(Pre.C5(:,indn),2)*size(Pre.C5(:,indn),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indn)))./(size(Pre.C0(:,indn),2)*size(Pre.C0(:,indn),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h5,'clim',[.45,.75] )
%h1.XLim = h1.XLim*0.9;
%h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_icm.fig'])
close all

%Next ...
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indn1),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn1),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indn1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indn1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indn1),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indn1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indn1),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn1),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn1),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn1),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indn1),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indn1),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indn1),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indn1),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn1)))./(size(Pre.C5(:,indn1),2)*size(Pre.C5(:,indn1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indn1)))./(size(Pre.C3(:,indn1),2)*size(Pre.C3(:,indn1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indn1)))./(size(Pre.C1(:,indn1),2)*size(Pre.C1(:,indn1),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h1,'clim',[.45,.75] )
set(h2,'clim',[.45,.75] )
set(h3,'clim',[.45,.75] )
set(h4,'clim',[.45,.75] )
set(h5,'clim',[.45,.75] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human_1'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_1.fig'])
close all


h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn1)))./(size(Pre.C5(:,indn1),2)*size(Pre.C5(:,indn1),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indn1)))./(size(Pre.C0(:,indn1),2)*size(Pre.C0(:,indn1),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h5,'clim',[.45,.75] )
%h1.XLim = h1.XLim*0.9;
%h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human_1_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_1_icm.fig'])
close all


%Next next next
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices', OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indn2),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn2),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indn2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indn2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indn2),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indn2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny
colormap(parula);

h4 = subplot(2,3,3);

f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indn2),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn2),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn2),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn2),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indn2),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indn2),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indn2),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indn2),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn2)))./(size(Pre.C5(:,indn2),2)*size(Pre.C5(:,indn2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indn2)))./(size(Pre.C3(:,indn2),2)*size(Pre.C3(:,indn2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indn2)))./(size(Pre.C1(:,indn2),2)*size(Pre.C1(:,indn2),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h1,'clim',[.45,.75] )
set(h2,'clim',[.45,.75] )
set(h3,'clim',[.45,.75] )
set(h4,'clim',[.45,.75] )
set(h5,'clim',[.45,.75] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human_2'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_2.fig'])
close all

h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn2)))./(size(Pre.C5(:,indn2),2)*size(Pre.C5(:,indn2),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indn2)))./(size(Pre.C0(:,indn2),2)*size(Pre.C0(:,indn2),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h5,'clim',[.45,.75] )
%h1.XLim = h1.XLim*0.9;
%h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human_2_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_2_icm.fig'])
close all

%Next
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
%view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn3),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn3),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices', OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indn3),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn3),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indn3),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indn3),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indn3),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indn3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indn3),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn3),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn3),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn3),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indn3),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indn3),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indn3),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indn3),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn3)))./(size(Pre.C5(:,indn3),2)*size(Pre.C5(:,indn3),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indn3)))./(size(Pre.C3(:,indn3),2)*size(Pre.C3(:,indn3),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indn3)))./(size(Pre.C1(:,indn3),2)*size(Pre.C1(:,indn3),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);
set(h1,'clim',[.45,.75] )
set(h2,'clim',[.45,.75] )
set(h3,'clim',[.45,.75] )
set(h4,'clim',[.45,.75] )
set(h5,'clim',[.45,.75] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;

print(['./Plots/Naive_Human_3'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_3.fig'])
close all


h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn3)))./(size(Pre.C5(:,indn3),2)*size(Pre.C5(:,indn3),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indn3)))./(size(Pre.C0(:,indn3),2)*size(Pre.C0(:,indn3),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h5,'clim',[.45,.75] )
%h1.XLim = h1.XLim*0.9;
%h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human_3_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_3_icm.fig'])
close all

%Next
h = figure(2)
h1 = subplot(2,3,5);
f2 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn4),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view(169.9208,90)
view([375.7415,-26.7022])
camlight('left')
material shiny
colormap(parula);

h2 = subplot(2,3,6);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn4),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn4),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',  OBJ2c.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn4),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([-70.8201,90])
camlight('left')
material shiny
colormap(parula);

h3 = subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1(:,indn4),2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2(:,indn4),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4(:,indn4),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5(:,indn4),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6(:,indn4),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7(:,indn4),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material shiny 
colormap(parula);

h4 = subplot(2,3,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1(:,indn4),2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn4),2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3(:,indn4),2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2(:,indn4),2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4(:,indn4),2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5(:,indn4),2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7(:,indn4),2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6(:,indn4),2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material shiny 
colormap(parula);

h5 = subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn4)))./(size(Pre.C5(:,indn4),2)*size(Pre.C5(:,indn4),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3(:,indn4)))./(size(Pre.C3(:,indn4),2)*size(Pre.C3(:,indn4),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1(:,indn4)))./(size(Pre.C1(:,indn4),2)*size(Pre.C1(:,indn4),1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material shiny 
colormap(parula);

set(h1,'clim',[.45,.75] )
set(h2,'clim',[.45,.75] )
set(h3,'clim',[.45,.75] )
set(h4,'clim',[.45,.75] )
set(h5,'clim',[.45,.75] )
h1.XLim = h1.XLim*0.9;
h1.YLim = h1.YLim*0.9;

print(['./Plots/Naive_Human_4'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_4.fig'])
close all

h = figure(2)
h5 = subplot(2,3,1);
f1 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C5(:,indn4)))./(size(Pre.C5(:,indn4),2)*size(Pre.C5(:,indn4),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',sum(sum(Pre.C0(:,indn4)))./(size(Pre.C0(:,indn4),2)*size(Pre.C0(:,indn4),1)*ones(length(O5.vertices),1)),'FaceColor','interp','LineStyle','none');
view([181.1552,-78.9270])
axis equal
axis off
camlight('left')
material shiny
colormap(parula);
set(h5,'clim',[.45,.75] )
%h1.XLim = h1.XLim*0.9;
%h1.YLim = h1.YLim*0.9;
print(['./Plots/Naive_Human_4_icm'],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Naive_Human_4_icm.fig'])
close all


