function h = PlotEmbryoCS7(varargin);

%OBJ1,figno,tissues,XYZ
if length(varargin)==2
    OBJ3 = varargin{1};
    tissues = varargin{2};
    figno = 1;
    XYZ = [];    
    Y = [];
    Yp = [];
    XYZp = [];
elseif length(varargin)==3    
    OBJ3 = varargin{1};
    tissues = varargin{2};
    figno = varargin{3};
    XYZ = [];    
    Y = [];
    Yp = [];
    XYZp = [];
elseif length(varargin)==4
    OBJ3 = varargin{1};
    tissues = varargin{2};
    figno = varargin{3};
    Output = varargin{4};
    XYZ = Output.cleanX ; 
    XYZp = [];
    Y = [];
elseif length(varargin)==5
    OBJ3 = varargin{1};
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
    tissues = {'EmDisc','Am','ExMes','Tb','SYS','VE','Stalk'};
end

if isempty(XYZp)==0
   opac = 0.1;
else
    opac = 1;
end

h = figure(figno);

if sum(strcmp(tissues,'Am'))~=0
f1 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices', OBJ3.vertices,'FaceColor',[26,8,115]/255,'LineStyle','none','FaceAlpha',opac);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'Am_CS7')==1 );
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'Am_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end
end

if sum(strcmp(tissues,'EmDisc'))~=0
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[2,51,191]/255,'LineStyle','none','FaceAlpha',opac);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'EmDisc_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'EmDisc_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end


if sum(strcmp(tissues,'Stalk'))~=0
f8 = patch('Faces',OBJ3.objects(30).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[96, 56, 19]/255,'LineStyle','none','FaceAlpha',opac);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'Stalk_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'Stalk_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end


if sum(strcmp(tissues,'VE'))~=0
f4 = patch('Faces',OBJ3.objects(12).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[191,60,4]/255,'LineStyle','none','FaceAlpha',opac);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'VE_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'VE_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end
end


if sum(strcmp(tissues,'SYS'))~=0
f5 = patch('Faces',OBJ3.objects(8).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[191,113,4]/255,'LineStyle','none','FaceAlpha',opac);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'SYS_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'SYS_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end




if sum(strcmp(tissues,'Tb'))~=0
f6 = patch('Faces',OBJ3.objects(4).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[113,8,166]/255,'LineStyle','none','FaceAlpha',opac);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'Tb_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'Tb_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end

end

if sum(strcmp(tissues,'ExMes'))~=0
f7 = patch('Faces',OBJ3.objects(34).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[150,119,0]/255,'LineStyle','none','FaceAlpha',opac);
hold on
    if isempty(XYZp)==0 & isempty(Yp)==1
        indT = find(strcmp(Output.cleanAnotaton,'ExMes_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
    elseif isempty(XYZp)==0 & isempty(Yp)==0
        indT = find(strcmp(Output.cleanAnotaton,'ExMes_CS7')==1);
        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')
    end
end

axis equal
axis off
camlight('right')
material dull 
colormap(parula);

output = 1;