function h = PlotEmbryoCS5GP(varargin)

if length(varargin)==4
    Output = varargin{1};
    OBJ1 = varargin{2};
    tissues = varargin{3};
    fig = varargin{4};
else
    Output = varargin{1};
    OBJ1 = varargin{2};
    tissues = varargin{3};
    fig = varargin{4};
    subplotp = varargin{5};
end

h = figure(fig)
if length(varargin)==5
    h1 = subplot(subplotp(1),subplotp(2),[subplotp(3:end)]);
end

if strcmp(tissues,'all')==1
    tissues = {'EmDisc','Am','ExMes','Tb','SYS','VE'}
end


if sum(strcmp(tissues,'Am'))~=0
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceVertexCData',Output.m1,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'EmDisc'))~=0
f1 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',  OBJ1.vertices,'FaceVertexCData',Output.m2,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'VE'))~=0
f2 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',  OBJ1.vertices,'FaceVertexCData',Output.m4,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'SYS'))~=0
f3 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',  OBJ1.vertices,'FaceVertexCData',Output.m5,'FaceColor','interp','LineStyle','none');
end


if sum(strcmp(tissues,'ExMes'))~=0
f4 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',  OBJ1.vertices,'FaceVertexCData',Output.m6,'FaceColor','interp','LineStyle','none');
end

if sum(strcmp(tissues,'Tb'))~=0
f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',  OBJ1.vertices,'FaceVertexCData',Output.m7,'FaceColor','interp','LineStyle','none');
end

axis equal
axis off
material dull 
colormap(parula);
%view([-56.3540,-8.2808])
%camlight('right')



output = 1;