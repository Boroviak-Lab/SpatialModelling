function [OBJout,a,b,OutputO] = transformCS62(varargin)

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

OBJT2 = transformOBJ(OBJ1,pi/4,[0,0,1]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(-OutputO.cleanX,pi/4,[0,0,1]);
    OutputO.cleanX = Quaternion3(pi/4,[0,0,1],OutputO.cleanX);

end

OBJT3 = transformOBJ(OBJT2,'flipX');
if length(varargin)==3
OutputO.cleanX(:,1) = -OutputO.cleanX(:,1);
end

OBJout = transformOBJ(OBJT3,'flipZ');
if length(varargin)==3
OutputO.cleanX(:,3) = -OutputO.cleanX(:,3);
end

a = 20.5492;
b = -60.9143;
%Transcformation done

else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation to get CS6 complete embryo view %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transform the EmDisc/VE views
%h = PlotEmbryoCS6GP(Output,OBJ1,{'EmDisc','Stalk','VE'},2);


OBJT2 = transformOBJ(OBJ1,pi/4,[0,0,1]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(-OutputO.cleanX,pi/4,[0,0,1]);
    OutputO.cleanX = Quaternion3(pi/4,[0,0,1],OutputO.cleanX);


end

OBJout = OBJT2;

%OBJout = transformOBJ(OBJT2,'flipX');
%if length(varargin)==2
%OutputO.cleanX(:,1) = -OutputO.cleanX(:,1);
%end
%h = PlotEmbryoCS6GP(Output,OBJT3,{'EmDisc','Stalk','VE'},4);
%view([0,90])

a = -209.4689;
b= -77.0435;

%a = 38.8691;
%b = -82.5458
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end