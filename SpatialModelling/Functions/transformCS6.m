function [OBJout,a,b,OutputO] = transformCS6(varargin)

if length(varargin)==2
    OBJ1 = varargin{1};
    mode = varargin{2};
    OutputO = [];
else
    OBJ1 = varargin{1};
    OutputO = varargin{2};
    mode = varargin{3};
end

if strcmp(mode,'all')==1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation to get CS6 complete embryo view %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa = linspace(0,-2*pi,100);

i = 11; %1:100
OBJT2 = transformOBJ(OBJ1,aa(i),[1,0,0]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(OutputO.cleanX,aa(i),[1,0,0]);
    OutputO.cleanX = Quaternion3(aa(i),[1,0,0],OutputO.cleanX);
end
%h = PlotEmbryoCS6GP(Output,OBJT2,{'all'},2);
%view([183.0698,90])
%pause
%clf
%end
%camlight('right')
%camlight('left')
%i = 11;
%view([162.4898,83.3904])
%view([177.9950,80.9869])

OBJT3 = transformOBJ(OBJT2,'flipX');
if length(varargin)==3
OutputO.cleanX(:,1) = -OutputO.cleanX(:,1);
end

%h = PlotEmbryoCS6GP(Output,OBJT3,{'EmDisc','Am','Tb','Stalk','VE','SYS'},3);
%view([177.9950,80.9869])
%camlight('left')
%view([164.8625,74.9781])
%view([131.1068,61.7588])


aa = linspace(0,-2*pi,100);
i = 5; %1:100
OBJout = transformOBJ(OBJT3,aa(i),[1,0,0]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(OutputO.cleanX,aa(i),[1,0,0]);
    OutputO.cleanX = Quaternion3(aa(i),[1,0,0],OutputO.cleanX);
end
%h = PlotEmbryoCS6GP(Output,OBJT4,{'all'},4);
%view([131.1068,61.7588])
%camlight('left')
%pause
%clf
%end

%h = PlotEmbryoCS6GP(Output,OBJT4,{'EmDisc','Am','Tb','Stalk','VE','SYS'},4);
%view([131.1068,61.7588])
%camlight('left')
a = 135.4947; 
b = 65.3641;
%Transcformation done

else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation to get CS6 complete embryo view %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transform the EmDisc/VE views
%h = PlotEmbryoCS6GP(Output,OBJ1,{'EmDisc','Stalk','VE'},2);

aa = linspace(0,2*pi,100);
i = 15; %1:100
OBJT1 = transformOBJ(OBJ1,aa(i),[0,1,0]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(OutputO.cleanX,aa(i),[0,1,0]);

    OutputO.cleanX = Quaternion3(aa(i),[0,1,0],OutputO.cleanX);

end
%h = PlotEmbryoCS6GP(Output,OBJT1,{'EmDisc','Stalk','VE'},2);
%view([0,90])
%pause
%clf
%end

i = 15%1:100

OBJT2 = transformOBJ(OBJT1,aa(i),[1,0,0]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(OutputO.cleanX,aa(i),[1,0,0]);
    OutputO.cleanX = Quaternion3(aa(i),[1,0,0],OutputO.cleanX);
end
%h = PlotEmbryoCS6GP(Output,OBJT2,{'EmDisc','Stalk','VE'},3);
%view([0,90])
%pause%
%clf
%end

OBJout = transformOBJ(OBJT2,'flipY');
if length(varargin)==3
OutputO.cleanX(:,2) = -OutputO.cleanX(:,2);
end
%h = PlotEmbryoCS6GP(Output,OBJT3,{'EmDisc','Stalk','VE'},4);
%view([0,90])
a = -113.7252;
b = -21.0065;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end