addpath(genpath('Functions'))

%Load CS7
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,ShotsCS5] = LoadShots('CS5');
%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');
%Load CS7
[OBJ3,section] = LoadCS7('3D');
[D,Locations,XYZ,CellType,ShotsCS7] = LoadShots('CS7');

%Load in the IDs for the marmoset correlation matrices (these were
%generated in R following loading in of datasets and alignments with human
%cells)
LOCMarm = importdata('../Data/SpatialData/PrimedNiave/EmbryoCorrHeatmaps/Cor_Loc1.csv')
IDMarm = importdata('../DataSpatialData/PrimedNiave/EmbryoCorrHeatmaps/Cor_ID1.csv')

%Correlations one for each batch of cells (3 conventional, 2 niave)
C = importdata('../DataSpatialData/PrimedNiave/Cor_Conv1.csv');

R5 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Tb_CS3')==1),:);
R3 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Hyp_CS3')==1),:);
R1 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Epi_CS3')==1),:);
R0 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'ICM_CS3')==1),:);

%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);
[OutputCS7] = loadCS7ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS7);

[OutputCS5] = MarmosetGP_CS5Corr(OutputCS5);
[OutputCS6] = MarmosetGP_CS6Corr(OutputCS6);
[OutputCS7] = MarmosetGP_CS7Corr(OutputCS7);

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

[OutputCS5] = MarmosetGPInfer_CS5Corr(OutputCS5,1,OBJ1);
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,1,OBJ2);
%[OutputCS7] = MarmosetGPInfer_CS7Corr(OutputCS7,1,OBJ3);

m1_CS5 = zeros(length(OutputCS5.m1),1);
m2_CS5 = zeros(length(OutputCS5.m2),1);
m3_CS5 = zeros(length(OutputCS5.m3),1);
m4_CS5 = zeros(length(OutputCS5.m4),1);
m5_CS5 = zeros(length(OutputCS5.m5),1);
m6_CS5 = zeros(length(OutputCS5.m6),1);
m7_CS5 = zeros(length(OutputCS5.m7),1);
m1_CS6 = zeros(length(OutputCS6.m1),1);
m2_CS6 = zeros(length(OutputCS6.m2),1);
m3_CS6 = zeros(length(OutputCS6.m3),1);
m4_CS6 = zeros(length(OutputCS6.m4),1);
m5_CS6 = zeros(length(OutputCS6.m5),1);
m6_CS6 = zeros(length(OutputCS6.m6),1);
m7_CS6 = zeros(length(OutputCS6.m7),1);
m8_CS6 = zeros(length(OutputCS6.m8),1);
%m1_CS7 = zeros(length(OutputCS7.m1),1);
%m2_CS7 = zeros(length(OutputCS7.m2),1);
%m3_CS7 = zeros(length(OutputCS7.m3),1)';
%m5_CS7 = zeros(length(OutputCS7.m5),1)';
%m6_CS7 = zeros(length(OutputCS7.m6),1)';
%m7_CS7 = zeros(length(OutputCS7.m7),1)';

IDIV = importdata('../DataSpatialData/PrimedNiave/Cor_uIDC1.csv')

for i = 1 :size(C.data,2)

[OutputCS5] = MarmosetGPInfer_CS5Corr(OutputCS5,i,OBJ1);
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,i,OBJ2);
%[OutputCS7] = MarmosetGPInfer_CS7Corr(OutputCS7,i,OBJ3);

m1_CS5 = [m1_CS5 , OutputCS5.m1];
m2_CS5 = [m2_CS5 , OutputCS5.m2];
m3_CS5 = [m3_CS5 , OutputCS5.m3];
m4_CS5 = [m4_CS5 , OutputCS5.m4];
m5_CS5 = [m5_CS5 , OutputCS5.m5];
m6_CS5 = [m6_CS5 , OutputCS5.m6];
m7_CS5 = [m7_CS5 , OutputCS5.m7];
m1_CS6 = [m1_CS6 , OutputCS6.m1];
m2_CS6 = [m2_CS6 , OutputCS6.m2];
m3_CS6 = [m3_CS6 , OutputCS6.m3];
m4_CS6 = [m4_CS6 , OutputCS6.m4];
m5_CS6 = [m5_CS6 , OutputCS6.m5];
m6_CS6 = [m6_CS6 , OutputCS6.m6];
m7_CS6 = [m7_CS6 , OutputCS6.m7];
m8_CS6 = [m8_CS6 , OutputCS6.m8];
% m1_CS7 = [m1_CS ; OutputCS7.m1];
% m2_CS7 = [m2_CS7 ; OutputCS7.m2];
% m3_CS7 = [m3_CS7 ; OutputCS7.m3];
% m5_CS7 = [m5_CS7 ; OutputCS7.m5];
% m6_CS7 = [m6_CS7 ; OutputCS7.m6];
% m7_CS7 = [m7_CS7 ; OutputCS7.m7];

h = figure(2)
subplot(2,3,5);
%f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices',  OBJ1c.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1f.objects(20).data.vertices,'Vertices',  OBJ1f.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([26.5143,-23.6498])
%view([c1,d1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )

subplot(2,3,6);
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
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
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,3);

f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R5(:,i)))./(size(R5,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R3(:,i)))./(size(R3,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R1(:,i)))./(size(R1,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
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

print(['./Plots/' filem{1}],'-dpdf','-r1000');
%savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS6_3D_' gene2test '.fig'])
clf
end

%Correlations for batch 2
C = importdata('../DataSpatialData/PrimedNiave/Cor_Conv2.csv');
IDIV = importdata('../DataSpatialData/PrimedNiave/Cor_uIDC2.csv')

R5_2 = [C.data(find(strcmp(IDMarm.textdata(2:end,2),'Tb_CS3')==1),:)];
R3_2 = [C.data(find(strcmp(IDMarm.textdata(2:end,2),'Hyp_CS3')==1),:)];
R1_2 = [C.data(find(strcmp(IDMarm.textdata(2:end,2),'Epi_CS3')==1),:)];


%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);
[OutputCS7] = loadCS7ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS7);

[OutputCS5] = MarmosetGP_CS5Corr(OutputCS5);
[OutputCS6] = MarmosetGP_CS6Corr(OutputCS6);
[OutputCS7] = MarmosetGP_CS7Corr(OutputCS7);

for i = 1 :size(C.data,2)

[OutputCS5] = MarmosetGPInfer_CS5Corr(OutputCS5,i,OBJ1);
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,i,OBJ2);
[OutputCS7] = MarmosetGPInfer_CS7Corr(OutputCS7,i,OBJ3);

m1_CS5 = [m1_CS5 , OutputCS5.m1];
m2_CS5 = [m2_CS5 , OutputCS5.m2];
m3_CS5 = [m3_CS5 , OutputCS5.m3];
m4_CS5 = [m4_CS5 , OutputCS5.m4];
m5_CS5 = [m5_CS5 , OutputCS5.m5];
m6_CS5 = [m6_CS5 , OutputCS5.m6];
m7_CS5 = [m7_CS5 , OutputCS5.m7];
m1_CS6 = [m1_CS6 , OutputCS6.m1];
m2_CS6 = [m2_CS6 , OutputCS6.m2];
m3_CS6 = [m3_CS6 , OutputCS6.m3];
m4_CS6 = [m4_CS6 , OutputCS6.m4];
m5_CS6 = [m5_CS6 , OutputCS6.m5];
m6_CS6 = [m6_CS6 , OutputCS6.m6];
m7_CS6 = [m7_CS6 , OutputCS6.m7];
m8_CS6 = [m8_CS6 , OutputCS6.m8];

%m1_CS7 = [m1_CS7 + OutputCS7.m1;
%m2_CS7 = [m2_CS7 + OutputCS7.m2;
%m3_CS7 = m3_CS7 + OutputCS7.m3;
%m5_CS7 = m5_CS7 + OutputCS7.m5;
%m6_CS7 = m6_CS7 + OutputCS7.m6;
%m7_CS7 = m7_CS7 + OutputCS7.m7;

h = figure(2)
subplot(2,3,5);
%f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices',  OBJ1c.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1f.objects(20).data.vertices,'Vertices',  OBJ1f.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([26.5143,-23.6498])
%view([c1,d1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )

subplot(2,3,6);
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
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
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,3);

f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R5_2(:,i)))./(size(R5_2,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R3_2(:,i)))./(size(R3_2,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R1_2(:,i)))./(size(R1_2,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
view([164.2550,-40.2050])
axis equal
axis off
%view([a2,b2])
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
filem = ['Conv2_' IDIV.textdata{i+1,2} '_' num2str(i) '_EmDAm.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/' filem{1}],'-dpdf','-r1000');
clf
end




%Correlations
C = importdata('../DataSpatialData/PrimedNiave/Cor_Conv3.csv');
IDIV = importdata('../DataSpatialData/PrimedNiave/Cor_uIDC3.csv')

R5_3 = [C.data(find(strcmp(IDMarm.textdata(2:end,2),'Tb_CS3')==1),:)];
R3_3 = [C.data(find(strcmp(IDMarm.textdata(2:end,2),'Hyp_CS3')==1),:)];
R1_3 = [C.data(find(strcmp(IDMarm.textdata(2:end,2),'Epi_CS3')==1),:)];

%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);
[OutputCS7] = loadCS7ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS7);

[OutputCS5] = MarmosetGP_CS5Corr(OutputCS5);
[OutputCS6] = MarmosetGP_CS6Corr(OutputCS6);
[OutputCS7] = MarmosetGP_CS7Corr(OutputCS7);


for i = 1 :size(C.data,2)

[OutputCS5] = MarmosetGPInfer_CS5Corr(OutputCS5,i,OBJ1);
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,i,OBJ2);
[OutputCS7] = MarmosetGPInfer_CS7Corr(OutputCS7,i,OBJ3);

m1_CS5 = [m1_CS5 , OutputCS5.m1];
m2_CS5 = [m2_CS5 , OutputCS5.m2];
m3_CS5 = [m3_CS5 , OutputCS5.m3];
m4_CS5 = [m4_CS5 , OutputCS5.m4];
m5_CS5 = [m5_CS5 , OutputCS5.m5];
m6_CS5 = [m6_CS5 , OutputCS5.m6];
m7_CS5 = [m7_CS5 , OutputCS5.m7];

m1_CS6 = [m1_CS6 , OutputCS6.m1];
m2_CS6 = [m2_CS6 , OutputCS6.m2];
m3_CS6 = [m3_CS6 , OutputCS6.m3];
m4_CS6 = [m4_CS6 , OutputCS6.m4];
m5_CS6 = [m5_CS6 , OutputCS6.m5];
m6_CS6 = [m6_CS6 , OutputCS6.m6];
m7_CS6 = [m7_CS6 , OutputCS6.m7];
m8_CS6 = [m8_CS6 , OutputCS6.m8];

%m1_CS7 = m1_CS7 + OutputCS7.m1;
%m2_CS7 = m2_CS7 + OutputCS7.m2;
%m3_CS7 = m3_CS7 + OutputCS7.m3;
%m5_CS7 = m5_CS7 + OutputCS7.m5;
%m6_CS7 = m6_CS7 + OutputCS7.m6;
%m7_CS7 = m7_CS7 + OutputCS7.m7;

h = figure(2)
subplot(2,3,5);
%f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices',  OBJ1c.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1f.objects(20).data.vertices,'Vertices',  OBJ1f.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([26.5143,-23.6498])
%view([c1,d1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )

subplot(2,3,6);
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
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
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,3);

f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R5_3(:,i)))./(size(R5_3,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R3_3(:,i)))./(size(R3_3,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R1_3(:,i)))./(size(R1_3,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
view([164.2550,-40.2050])
axis equal
axis off
%view([a2,b2])
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
filem = ['Conv3_' IDIV.textdata{i+1,2} '_' num2str(i) '_EmDAm.pdf']
%filem = ['Conv3_' num2str(i) '_EmDAm.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/' filem{1}],'-dpdf','-r1000');
clf
end

%We can save these all, but will be large files
MappingCS5.m1 = m1_CS5(:,2:end);
MappingCS5.m2 = m2_CS5(:,2:end);
MappingCS5.m3 = m3_CS5(:,2:end);
MappingCS5.m4 = m4_CS5(:,2:end);
MappingCS5.m5 = m5_CS5(:,2:end);
MappingCS5.m6 = m6_CS5(:,2:end);
MappingCS5.m7 = m7_CS5(:,2:end);
%save(['./Plots/CS5Mapping.mat'],'MappingCS5')

MappingCS6.m1 = m1_CS6(:,2:end);
MappingCS6.m2 = m2_CS6(:,2:end);
MappingCS6.m3 = m3_CS6(:,2:end);
MappingCS6.m4 = m4_CS6(:,2:end);
MappingCS6.m5 = m5_CS6(:,2:end);
MappingCS6.m6 = m6_CS6(:,2:end);
MappingCS6.m7 = m7_CS6(:,2:end);
MappingCS6.m8 = m8_CS6(:,2:end);
%save(['./Plots/CS6Mapping.mat'],'MappingCS6')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now do niave correlations
C = importdata('../DataSpatialData/PrimedNiave/Cor_Naive1.csv');
IDIV = importdata('../DataSpatialData/PrimedNiave/Cor_uIDN1.csv')

R5_5 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Tb_CS3')==1),:);
R3_5 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Hyp_CS3')==1),:);
R1_5 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Epi_CS3')==1),:);

%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);
%[OutputCS7] = loadCS7ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS7);

[OutputCS5] = MarmosetGP_CS5Corr(OutputCS5);
[OutputCS6] = MarmosetGP_CS6Corr(OutputCS6);
%[OutputCS7] = MarmosetGP_CS7Corr(OutputCS7);

[OutputCS5] = MarmosetGPInfer_CS5Corr(OutputCS5,1,OBJ1);
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,1,OBJ2);
%[OutputCS7] = MarmosetGPInfer_CS7Corr(OutputCS7,1,OBJ3);

m1_CS5 = zeros(length(OutputCS5.m1),1);
m2_CS5 = zeros(length(OutputCS5.m2),1);
m3_CS5 = zeros(length(OutputCS5.m3),1);
m4_CS5 = zeros(length(OutputCS5.m4),1);
m5_CS5 = zeros(length(OutputCS5.m5),1);
m6_CS5 = zeros(length(OutputCS5.m6),1);
m7_CS5 = zeros(length(OutputCS5.m7),1);
m1_CS6 = zeros(length(OutputCS6.m1),1);
m2_CS6 = zeros(length(OutputCS6.m2),1);
m3_CS6 = zeros(length(OutputCS6.m3),1);
m4_CS6 = zeros(length(OutputCS6.m4),1);
m5_CS6 = zeros(length(OutputCS6.m5),1);
m6_CS6 = zeros(length(OutputCS6.m6),1);
m7_CS6 = zeros(length(OutputCS6.m7),1);

for i = 1 :size(C.data,2)


[OutputCS5] = MarmosetGPInfer_CS5Corr(OutputCS5,i,OBJ1);
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,i,OBJ2);
[OutputCS7] = MarmosetGPInfer_CS7Corr(OutputCS7,i,OBJ3);

m1_CS5 = [m1_CS5 , OutputCS5.m1];
m2_CS5 = [m2_CS5 , OutputCS5.m2];
m3_CS5 = [m3_CS5 , OutputCS5.m3];
m4_CS5 = [m4_CS5 , OutputCS5.m4];
m5_CS5 = [m5_CS5 , OutputCS5.m5];
m6_CS5 = [m6_CS5 , OutputCS5.m6];
m7_CS5 = [m7_CS5 , OutputCS5.m7];
m1_CS6 = [m1_CS6 , OutputCS6.m1];
m2_CS6 = [m2_CS6 , OutputCS6.m2];
m3_CS6 = [m3_CS6 , OutputCS6.m3];
m4_CS6 = [m4_CS6 , OutputCS6.m4];
m5_CS6 = [m5_CS6 , OutputCS6.m5];
m6_CS6 = [m6_CS6 , OutputCS6.m6];
m7_CS6 = [m7_CS6 , OutputCS6.m7];
m8_CS6 = [m8_CS6 , OutputCS6.m8];
%m2_CS7 = m2_CS7 + OutputCS7.m2;
%m3_CS7 = m3_CS7 + OutputCS7.m3;
%m5_CS7 = m5_CS7 + OutputCS7.m5;
%m6_CS7 = m6_CS7 + OutputCS7.m6;
%m7_CS7 = m7_CS7 + OutputCS7.m7;

h = figure(2)
subplot(2,3,5);
%f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices',  OBJ1c.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1f.objects(20).data.vertices,'Vertices',  OBJ1f.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([26.5143,-23.6498])
%view([c1,d1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )

subplot(2,3,6);
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
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
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,3);

f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R5_5(:,i)))./(size(R5_5,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R3_5(:,i)))./(size(R3_5,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R1_5(:,i)))./(size(R1_5,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
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

filem = ['Naive1_' num2str(i) '_EmDAm.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/' filem{1}],'-dpdf','-r1000');
clf
end

%Correlations
C = importdata('../DataSpatialData/PrimedNiave/Cor_Naive2.csv');
IDIV = importdata('../DataSpatialData/PrimedNiave/Cor_uIDN2.csv')

R5_6 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Tb_CS3')==1),:);
R3_6 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Hyp_CS3')==1),:);
R1_6 = C.data(find(strcmp(IDMarm.textdata(2:end,2),'Epi_CS3')==1),:);

%Proecss CS6 data
[OutputCS5] = loadCS5ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS5);
[OutputCS6] = loadCS6ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS6);
%[OutputCS7] = loadCS7ScaffoldCorr(C,LOCMarm,IDMarm,ShotsCS7);

[OutputCS5] = MarmosetGP_CS5Corr(OutputCS5);
[OutputCS6] = MarmosetGP_CS6Corr(OutputCS6);
%[OutputCS7] = MarmosetGP_CS7Corr(OutputCS7);


for i = 1 :size(C.data,2)

[OutputCS5] = MarmosetGPInfer_CS5Corr(OutputCS5,i,OBJ1);
[OutputCS6] = MarmosetGPInfer_CS6Corr(OutputCS6,i,OBJ2);
%[OutputCS7] = MarmosetGPInfer_CS7Corr(OutputCS7,i,OBJ3);

m1_CS5 = [m1_CS5 , OutputCS5.m1];
m2_CS5 = [m2_CS5 , OutputCS5.m2];
m3_CS5 = [m3_CS5 , OutputCS5.m3];
m4_CS5 = [m4_CS5 , OutputCS5.m4];
m5_CS5 = [m5_CS5 , OutputCS5.m5];
m6_CS5 = [m6_CS5 , OutputCS5.m6];
m7_CS5 = [m7_CS5 , OutputCS5.m7];
m1_CS6 = [m1_CS6 , OutputCS6.m1];
m2_CS6 = [m2_CS6 , OutputCS6.m2];
m3_CS6 = [m3_CS6 , OutputCS6.m3];
m4_CS6 = [m4_CS6 , OutputCS6.m4];
m5_CS6 = [m5_CS6 , OutputCS6.m5];
m6_CS6 = [m6_CS6 , OutputCS6.m6];
m7_CS6 = [m7_CS6 , OutputCS6.m7];
m8_CS6 = [m8_CS6 , OutputCS6.m8];

%m1_CS7 = m1_CS7 + OutputCS7.m1;
%m2_CS7 = m2_CS7 + OutputCS7.m2;
%m3_CS7 = m3_CS7 + OutputCS7.m3;
%m5_CS7 = m5_CS7 + OutputCS7.m5;
%m6_CS7 = m6_CS7 + OutputCS7.m6;
%m7_CS7 = m7_CS7 + OutputCS7.m7;

h = figure(2)
subplot(2,3,5);
%f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices',  OBJ1c.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1f.objects(20).data.vertices,'Vertices',  OBJ1f.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([26.5143,-23.6498])
%view([c1,d1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )

subplot(2,3,6);
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
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
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1,b1])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,3);

f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
f8 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',OBJ2b.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a2,b2])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.5,.8] )


subplot(2,3,1);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R5_6(:,i)))./(size(R5_6,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R3_6(:,i)))./(size(R3_6,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum((R1_6(:,i)))./(size(R1_6,1)*ones(length(OBJ0b.vertices),1)),'FaceColor','interp','LineStyle','none');
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

filem = ['Naive2_' num2str(i) '_EmDAm.pdf']
filem = strrep(filem,"/","_")
print(['./Plots/' filem{1}],'-dpdf','-r1000');
clf
end

%Again, these will be big files so commented out
NPre.C5 = [R5_5,R5_6];
NPre.C3 = [R3_5,R3_6];
NPre.C1 = [R1_5,R1_6];
%save(['./Plots/CS5NPre.mat'],'NPre')

NMappingCS5.m1 = m1_CS5(:,2:end);
NMappingCS5.m2 = m2_CS5(:,2:end);
NMappingCS5.m3 = m3_CS5(:,2:end);
NMappingCS5.m4 = m4_CS5(:,2:end);
NMappingCS5.m5 = m5_CS5(:,2:end);
NMappingCS5.m6 = m6_CS5(:,2:end);
NMappingCS5.m7 = m7_CS5(:,2:end);
%save(['./Plots/CS5NMapping.mat'],'NMappingCS5')

NMappingCS6.m1 = m1_CS6(:,2:end);
NMappingCS6.m2 = m2_CS6(:,2:end);
NMappingCS6.m3 = m3_CS6(:,2:end);
NMappingCS6.m4 = m4_CS6(:,2:end);
NMappingCS6.m5 = m5_CS6(:,2:end);
NMappingCS6.m6 = m6_CS6(:,2:end);
NMappingCS6.m7 = m7_CS6(:,2:end);
NMappingCS6.m8 = m8_CS6(:,2:end);
%save(['./Plots/CS6NMapping.mat'],'NMappingCS6')