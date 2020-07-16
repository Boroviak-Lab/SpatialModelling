addpath(genpath('Functions'))

D3 = importdata('Data/Keycorrect_CPAll2.csv')
Dexp = importdata('Data/NormData.csv');

[DexpCS5,Arc,Xtest,OutputCS5]  = loadCS5Section_3_surf(D3,Dexp);
[DexpCS6,Arc,Xtest,OutputCS6]  = loadCS6Section_2_surf(D3,Dexp);
[DexpCS7,Arc,Xtest,OutputCS7]  = loadCS7Section_2_surf(D3,Dexp);

genelist1{1} = 'TFAP2A';
genelist1{2} = 'GATA3';
genelist1{3} = 'MIXL1';
genelist1{4} = 'CDX2';

for i = 1:length(genelist1)
    try

[Output1] = Marmoset3D_CS5_surface(DexpCS5,OutputCS5,'Base2',genelist1{i});
[Output2] = Marmoset3D_CS6_surface(DexpCS6,OutputCS6,'Base2',genelist1{i});
[Output3] = Marmoset3D_CS7_surface(DexpCS7,OutputCS7,'Base2',genelist1{i});

maxV = max([Output1.m_0;Output1.m_1;Output1.m_2;Output1.m_3;Output1.m_4;Output1.m_5;   Output2.m_0;Output2.m_1;Output2.m_2; Output2.m_3; Output2.m_4;Output2.m_5;Output2.m_6;    Output3.m_0;Output3.m_1;Output3.m_2;Output3.m_3;Output3.m_5]);
minV = min([Output1.m_0;Output1.m_1;Output1.m_2;Output1.m_3;Output1.m_4;Output1.m_5;   Output2.m_0;Output2.m_1;Output2.m_2; Output2.m_3; Output2.m_4;Output2.m_5;Output2.m_6;    Output3.m_0;Output3.m_1;Output3.m_2;Output3.m_3;Output3.m_5]);

h = figure(i)

OBJ1=importdata('Data/CS5_section_cs4.mat');
OBJ1.vertices = Quaternion3(-0.0349066,[0,0,1],OBJ1.vertices);

subplot(1,3,1);
title([genelist1{i}])
f1 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_1,'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_0,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_5,'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_4,'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_2,'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_3,'FaceColor','interp','LineStyle','none');
axis equal
axis off
set(gca,'clim',[minV maxV] )
plot(linspace(-200,200,1000),-1000*ones(1,1000),'k-','LineWidth',2)
set(gca,'clim',[minV maxV] )
set(gca,'fontsize', 12)
set(gca,'TickLength',[0 0])
set(gca,'linewidth',.7)

OBJ1=importdata('Data/CS6_section_cs2.mat');
OBJ1.vertices = Quaternion3(-0.95,[0,0,1],OBJ1.vertices);

subplot(1,3,2);
f11 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output2.m_1,'FaceColor','interp','LineStyle','none');
hold on
f21 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output2.m_0,'FaceColor','interp','LineStyle','none');
f31 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output2.m_3,'FaceColor','interp','LineStyle','none');
f41 = patch('Faces',OBJ1.objects(32).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output2.m_4,'FaceColor','interp','LineStyle','none');
f51 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output2.m_2,'FaceColor','interp','LineStyle','none');
f61 = patch('Faces',OBJ1.objects(28).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output2.m_5,'FaceColor','interp','LineStyle','none');
f71 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output2.m_6,'FaceColor','interp','LineStyle','none');
axis equal
axis off
set(gca,'clim',[minV maxV] )
plot(linspace(-200,200,1000),-1000*ones(1,1000),'k-','LineWidth',2)

OBJ1=importdata('Data/CS7_section_cs2.mat');
OBJ1.vertices = Quaternion3(-1.75,[0,0,1],OBJ1.vertices);

subplot(1,3,3);
f12 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output3.m_1,'FaceColor','interp','LineStyle','none');
hold on
f22 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output3.m_0,'FaceColor','interp','LineStyle','none');
f32 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output3.m_1,'FaceColor','interp','LineStyle','none');
f42 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output3.m_2,'FaceColor','interp','LineStyle','none');
f52 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output3.m_5,'FaceColor','interp','LineStyle','none');
axis equal
axis off
plot(linspace(-200,200,1000),-1000*ones(1,1000),'k-','LineWidth',2)
set(gca,'clim',[minV maxV] )
cb=colorbar;
cb.Position = cb.Position + [0.1 0.15 -0 -0.3]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)
set(gcf,'color','w');

%print(h,['~/Desktop/Sections/' genelist1{i} '_3.png'],'-dpng','-r2000'); 
%clf

    catch
        clf
    end

end
