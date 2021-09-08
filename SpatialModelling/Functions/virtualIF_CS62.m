function output = virtualIF_CS62(Output,OBJ2A,tissues,fig)

h = figure(fig)

if strcmp(tissues,'all')==1
    tissues = {'EmDisc','Am','ExMes','Tb','SYS','VE','Stalk'}
end

if sum(strcmp(tissues,'Am'))~=0
f1 = patch('Faces',OBJ2A.objects(4).data.vertices,'Vertices',  OBJ2A.vertices,'FaceVertexCData',Output.m1,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'EmDisc'))~=0
f2 = patch('Faces',OBJ2A.objects(8).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m2,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'Stalk'))~=0
f4 = patch('Faces',OBJ2A.objects(20).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m2,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'VE'))~=0
f5 = patch('Faces',OBJ2A.objects(12).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m4,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'SYS'))~=0
f6 = patch('Faces',OBJ2A.objects(16).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m5,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'Tb'))~=0
f7 = patch('Faces',OBJ2A.objects(24).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m5,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'ExMes'))~=0
f8 = patch('Faces',OBJ2A.objects(28).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m6,'FaceColor','interp','LineStyle','none');
end

axis equal
axis off
material dull 
colormap(parula);
%view([-56.3540,-8.2808])
camlight('right')



output = 1;