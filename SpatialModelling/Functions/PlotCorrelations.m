addpath(genpath('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/code/gpml-matlab-v3.6-2015-07-07'))

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
LOCMarm = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/EmbryoCorrHeatmaps/Cor_Loc1.csv')
IDMarm = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/EmbryoCorrHeatmaps/Cor_ID1.csv')
C = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/EmbryoCorrHeatmaps/Cor_Conv1.csv');

%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);
[OutputCS7] = loadCS7ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS7);

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');

OBJ1d = transformOBJ(OBJ1,pi/2,[0,1,0]);
OBJ1e = transformOBJ(OBJ1d,'flipX');
OBJ1f = transformOBJ(OBJ1e,'flipZ');

%Also load the blastocyst
OBJ0a=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/Blastocyst_halfcut_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ0b=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/Blastocyst_slice_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')

load(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS5Pre.mat'])
load(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS5Mapping.mat'])
meh=load(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS6Mapping.mat'])

h = figure(2)
subplot(2,3,5);
%f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices',  OBJ1c.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1f.objects(20).data.vertices,'Vertices',  OBJ1f.vertices,'FaceVertexCData',mean(MappingCS5.m2,2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([26.5143,-23.6498])
%view([c1,d1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )

subplot(2,3,6);
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2,2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m8,2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3,2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([130.5487,90])
%view([-100.6086,15.6542])
%view([132.7748,-12.8427])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )

subplot(2,3,2);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m1,2),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m2,2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m4,2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m5,2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m6,2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(MappingCS5.m7,2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,3);

f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m1,2),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m2,2),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m3,2),'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m8,2),'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m4,2),'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m5,2),'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m7,2),'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',mean(MappingCS6.m6,2),'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C5))./(size(Pre.C5,2)*size(Pre.C5,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C3))./(size(Pre.C3,2)*size(Pre.C3,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(Pre.C1))./(size(Pre.C1,2)*size(Pre.C1,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )
%cb=colorbar;
%cb.Position = cb.Position + [0.1 0.15 -0 -0.3]
%cy=get(cb,'YTick');
%set(cb,'YTick',[cy(1),cy(3),cy(end)])
%set(gca,'fontsize', 12)
%set(gcf,'color','w');
%axis equal
%axis off
%camlightcam

filem = ['Conv1_' IDIV.textdata{i+1,2} '_' num2str(i) '_EmDAm.pdf']
filem = strrep(filem,"/","_")

print(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/' filem{1}],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS6_3D_' gene2test '.fig'])
clf

