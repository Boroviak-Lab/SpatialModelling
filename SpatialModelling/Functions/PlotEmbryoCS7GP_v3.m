function h = PlotEmbryoCS7GP_v3(varargin)

if length(varargin)==4
    Output = varargin{1};
    OBJ2A = varargin{2};
    tissues = varargin{3};
    fig = varargin{4};
else
    Output = varargin{1};
    OBJ2A = varargin{2};
    tissues = varargin{3};
    fig = varargin{4};
    subplotp = varargin{5};
end

h = figure(fig)
if length(varargin)==5
    h2 = subplot(subplotp(1),subplotp(2),[subplotp(3:end)]);
end

if strcmp(tissues,'all')==1
    tissues = {'EmDisc','Am','ExMes','Tb','SYS','VE','Stalk'}
end

if sum(strcmp(tissues,'Am'))~=0
f1 = patch('Faces',OBJ2A.objects(16).data.vertices,'Vertices',  OBJ2A.vertices,'FaceVertexCData',Output.m1,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'EmDisc'))~=0
f2 = patch('Faces',OBJ2A.objects(20).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m2,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'Stalk'))~=0
f4 = patch('Faces',OBJ2A.objects(30).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m8,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'VE'))~=0
f5 = patch('Faces',OBJ2A.objects(12).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m4,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'SYS'))~=0
f6 = patch('Faces',OBJ2A.objects(8).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m5,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'Tb'))~=0
f7 = patch('Faces',OBJ2A.objects(4).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m7,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'ExMes'))~=0
f8 = patch('Faces',OBJ2A.objects(34).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m6,'FaceColor','interp','LineStyle','none');
end

axis equal
axis off
material dull 
colormap(parula);
%view([-56.3540,-8.2808])
camlight('right')



output = 1;