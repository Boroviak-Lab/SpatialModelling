function h = PlotVirtualIF(varargin);

%OBJ1,figno,tissues,XYZ
if length(varargin)==2
    OBJ1 = varargin{1};
    tissues = varargin{2};
    figno = 1;
    XYZ = [];    
    Y = [];
    XYZp = [];    
    Yp = [];
elseif length(varargin)==3    
    OBJ1 = varargin{1};
    tissues = varargin{2};
    figno = varargin{3};
    XYZ = [];    
    Y = [];
    XYZp = [];    
    Yp = [];    
elseif length(varargin)==4
    OBJ1 = varargin{1};
    tissues = varargin{2};
    figno = varargin{3};
    Output = varargin{4};
    XYZ = Output.cleanX ; 
    XYZp = [];
    Y = [];
    Yp = [];
elseif length(varargin)==5
    OBJ1 = varargin{1};
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
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',opac);
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
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',opac);
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
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',opac);
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
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',opac);
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
f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',opac);
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
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,200,0]/255,'LineStyle','none','FaceAlpha',opac);
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