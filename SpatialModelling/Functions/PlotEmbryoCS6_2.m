function h = PlotEmbryoCS6_2(varargin);

%OBJ1,figno,tissues,XYZ
if length(varargin)==2
    OBJ2B = varargin{1};
    tissues = varargin{2};
    figno = 1;
    XYZ = [];    
    Y = [];
    Yp = [];
    XYZp = [];
elseif length(varargin)==3    
    OBJ2B = varargin{1};
    tissues = varargin{2};
    figno = varargin{3};
    XYZ = [];    
    Y = [];
    Yp = [];
    XYZp = [];
elseif length(varargin)==4
    OBJ2B = varargin{1};
    tissues = varargin{2};
    figno = varargin{3};
    Output = varargin{4};
    XYZ = Output.cleanX ; 
    XYZp = [];
    Y = [];
elseif length(varargin)==5
    OBJ2B = varargin{1};
    tissues = varargin{2};
    figno = varargin{3};
    Output = varargin{4};
    XYZ = Output.cleanX ;  
    type = varargin{5};
    XYZp = 1;
    if strcmp(type,'Expression')==1
        Y = Output.Ytrain;
        Yp = 1;
    else
        Yp = [];
    end
else
    disp('Wrong number of input parameters: requires 2, 3, 4, or 5 inputs')
end

if strcmp(tissues,'all')==1
    tissues = {'EmDisc','Am','ExMes','Tb','SYS','VE'};
end

if isempty(XYZp)==0
   opac = 0.1;
else
    opac = 1;
end

h = figure(figno);

if sum(strcmp(tissues,'Am'))~=0
f1 = patch('Faces',OBJ2B.objects(4).data.vertices,'Vertices',  OBJ2B.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',1);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'Am_CS5')==1 | strcmp(Output.cleanAnotaton,'Am_CS5_PGC')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'Am_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end
end

if sum(strcmp(tissues,'EmDisc'))~=0
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceColor',[218,122,7]/255,'LineStyle','none','FaceAlpha',1);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end


if sum(strcmp(tissues,'Stalk'))~=0
f4 = patch('Faces',OBJ2B.objects(20).data.vertices,'Vertices',OBJ2B.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',1);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end

if sum(strcmp(tissues,'VE'))~=0
f5 = patch('Faces',OBJ2B.objects(12).data.vertices,'Vertices',OBJ2B.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',1);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'VE_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'VE_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end
end

if sum(strcmp(tissues,'SYS'))~=0
f6 = patch('Faces',OBJ2B.objects(16).data.vertices,'Vertices',OBJ2B.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',1);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'SYS_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'SYS_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end

if sum(strcmp(tissues,'Tb'))~=0
f7 = patch('Faces',OBJ2B.objects(24).data.vertices,'Vertices',OBJ2B.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',1);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'Tb_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'Tb_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end

if sum(strcmp(tissues,'ExMes'))~=0
f8 = patch('Faces',OBJ2B.objects(28).data.vertices,'Vertices',OBJ2B.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',1);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'ExMes_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'ExMes_CS5')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end
end

axis equal
axis off
camlight('right')
material dull 
colormap(parula);

output = 1;