addpath(genpath('Functions'))
%Load CS5
load('../Data/SpatialData/CS5_revisions_HIGHRES.mat','OBJ1')
load('../Data/SpatialData/Cross Sections/CS5/CS5_revisions_DIAGONAL_2.mat')

%Load CS6
load('../Data/SpatialData/Cross sections/CS6/CS6_revisions_CROSS_246.mat')
%load('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/CS6 Original Twin (EII)/05 Blender surfaces/CS6_embryomodel_260421_HIGHRES.mat')
[OBJ2A_2, sectcs6] = LoadCS6('246');

%CS62
load('../Data/SpatialData/CS6_EI_twin_30042_HIGHRES.mat')
[OBJ2B_3,sections] = LoadCS6_2('D2')

%Load CS7
load('../Data/SpatialData/200319_CS7_build_full_stalkcut_FORDYLAN.mat')
[OBJ3_2, sectcs6] = LoadCS7('C2');

[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,Shots);
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS6');
[OutputCS6] = loadCS6Scaffold(D,Locations,Shots);
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS62');
[OutputCS62] = loadCS6Scaffold2(D,Locations,Shots);
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS7');
[OutputCS7] = loadCS7Scaffold(D,Locations,Shots);

OBJ1_6b = transformOBJ(OBJ1_6,pi,[0,0,1]);
OBJ1_6b = transformOBJ(OBJ1_6b,'flipX');


OBJ2A_2b = transformOBJ(OBJ2A_2,0.81*pi,[0,0,1]);
OBJ2A_2b = transformOBJ(OBJ2A_2b,'flipX');

OBJ2B_3b = transformOBJ(OBJ2B_3,1.4*pi,[0,0,1]);

OBJ3_2b = transformOBJ(OBJ3_2,-0.55*pi,[0,0,1]);
OBJ3_2b = transformOBJ(OBJ3_2b,'flipX');

genelist = unique({'SOX2','T'});


%Plot the sections for the standard model
for i = 1:length(genelist)
    try
gene2test = genelist{i};
    
%CS5
[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,gene2test);
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1_6);

h1 = subplot(2,4,1);
%Amnion / 1 =Am/PGC
f1 = patch('Faces',OBJ1_6b.objects(4).data.vertices,'Vertices',  OBJ1_6b.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
%EmDisc / 2 = EmDisc
f1 = patch('Faces',OBJ1_6b.objects(20).data.vertices,'Vertices',  OBJ1_6b.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
%VE / 4 = VE
f2 = patch('Faces',OBJ1_6b.objects(8).data.vertices,'Vertices',  OBJ1_6b.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
%SYS / 5 = SYS
f3 = patch('Faces',OBJ1_6b.objects(12).data.vertices,'Vertices',  OBJ1_6b.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
%ExMes / 6 = ExMes
f4 = patch('Faces',OBJ1_6b.objects(16).data.vertices,'Vertices',  OBJ1_6b.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
%Tb / 7 = Tb
f5 = patch('Faces',OBJ1_6b.objects(24).data.vertices,'Vertices',  OBJ1_6b.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('right')
material dull 
colormap(parula);
title([gene2test])

%CS6
[OutputCS6] = MarmosetGP_CS6_v3(D,OutputCS6,gene2test);
[OutputCS6] = MarmosetGPInfer_CS6_v3(OutputCS6,OBJ2A_2);

%OBJ2A_2b = transformOBJ(OBJ2A_2,0.81*pi,[0,0,1]);

h2 = subplot(2,4,2);
%Am / 
f1 = patch('Faces',OBJ2A_2b.objects(20).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
%EmDisc
f2 = patch('Faces',OBJ2A_2b.objects(4).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
%Stalk
f3 = patch('Faces',OBJ2A_2b.objects(24).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m8,'FaceColor','interp','LineStyle','none');
%PGC
f4 = patch('Faces',OBJ2A_2b.objects(8).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
%VE
f5 = patch('Faces',OBJ2A_2b.objects(16).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
%SYS
f6 = patch('Faces',OBJ2A_2b.objects(12).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
%Tb
f7 = patch('Faces',OBJ2A_2b.objects(28).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
%ExMes
f8 = patch('Faces',OBJ2A_2b.objects(32).data.vertices,'Vertices',  OBJ2A_2b.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('right')
material dull 
colormap(parula);

[OutputCS62] = MarmosetGP_CS6(D,OutputCS62,gene2test);
[OutputCS62] = MarmosetGPInfer_CS6(OutputCS62,OBJ2B_3);    


h3 = subplot(2,4,4);
%Am
f1 = patch('Faces',OBJ2B_3b.objects(4).data.vertices,'Vertices',  OBJ2B_3b.vertices,'FaceVertexCData',OutputCS62.m1,'FaceColor','interp','LineStyle','none');
%EmDisc
f2 = patch('Faces',OBJ2B_3b.objects(8).data.vertices,'Vertices',  OBJ2B_3b.vertices,'FaceVertexCData',OutputCS62.m2,'FaceColor','interp','LineStyle','none');
%Stalk
f3 = patch('Faces',OBJ2B_3b.objects(20).data.vertices,'Vertices',  OBJ2B_3b.vertices,'FaceVertexCData',OutputCS62.m2,'FaceColor','interp','LineStyle','none');
%VE
f4 = patch('Faces',OBJ2B_3b.objects(12).data.vertices,'Vertices',  OBJ2B_3b.vertices,'FaceVertexCData',OutputCS62.m4,'FaceColor','interp','LineStyle','none');
%SYS
f5 = patch('Faces',OBJ2B_3b.objects(16).data.vertices,'Vertices',  OBJ2B_3b.vertices,'FaceVertexCData',OutputCS62.m5,'FaceColor','interp','LineStyle','none');
%Tb
f6 = patch('Faces',OBJ2B_3b.objects(24).data.vertices,'Vertices',  OBJ2B_3b.vertices,'FaceVertexCData',OutputCS62.m7,'FaceColor','interp','LineStyle','none');
%ExMes
f7 = patch('Faces',OBJ2B_3b.objects(28).data.vertices,'Vertices',  OBJ2B_3b.vertices,'FaceVertexCData',OutputCS62.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('right')
material dull 
colormap(parula);

%CS7
[OutputCS7] = MarmosetGP_CS7_v3(D,OutputCS7,gene2test);
[OutputCS7] = MarmosetGPInfer_CS7_v3(OutputCS7,OBJ3_2);

h4 = subplot(2,4,3);
%Am
f1 = patch('Faces',OBJ3_2b.objects(16).data.vertices,'Vertices',  OBJ3_2b.vertices,'FaceVertexCData',OutputCS7.m1,'FaceColor','interp','LineStyle','none');
%EmDisc
f2 = patch('Faces',OBJ3_2b.objects(20).data.vertices,'Vertices',  OBJ3_2b.vertices,'FaceVertexCData',OutputCS7.m2,'FaceColor','interp','LineStyle','none');
%Stalk
f3 = patch('Faces',OBJ3_2b.objects(30).data.vertices,'Vertices',  OBJ3_2b.vertices,'FaceVertexCData',OutputCS7.m8,'FaceColor','interp','LineStyle','none');
%SYS
f5 = patch('Faces',OBJ3_2b.objects(8).data.vertices,'Vertices',  OBJ3_2b.vertices,'FaceVertexCData',OutputCS7.m5,'FaceColor','interp','LineStyle','none');
%Tb
f6 = patch('Faces',OBJ3_2b.objects(4).data.vertices,'Vertices',  OBJ3_2b.vertices,'FaceVertexCData',OutputCS7.m7,'FaceColor','interp','LineStyle','none');
%ExMes
f7 = patch('Faces',OBJ3_2b.objects(34).data.vertices,'Vertices',  OBJ3_2b.vertices,'FaceVertexCData',OutputCS7.m6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('right')
material dull 
colormap(parula);

m1 = max(max([OutputCS5.cLim(2);OutputCS6.cLim(2);OutputCS7.cLim(2)]));
m2 = min(min([OutputCS5.cLim(1);OutputCS6.cLim(1);OutputCS7.cLim(1)]));


m1 = max([m1,0.07]);
m2 = max([m2,0]);

set(h1,'clim',[m2 m1] )
set(h2,'clim',[m2 m1] )
set(h3,'clim',[m2 m1] )
set(h4,'clim',[m2 m1] )

cb=colorbar;
cb.Position = cb.Position + [0.1 0.06 -0.0059 -0.0640]
cb.Color = 'k';
cb.FontSize = 12;
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
    try
set(cb,'YTick',[cy(1),cy(2),cy(end)])
    catch
set(cb,'YTick',[cy(1),cy(end)])        
    end
end


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
print('-dpng',['./Plots/' gene2test '_Section.png'],'-r1000')


clf
end
end