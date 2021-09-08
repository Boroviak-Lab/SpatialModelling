function [OBJout,a,b,OutputO] = transformCS7(varargin)

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
%h = PlotEmbryoCS7GP(Output,OBJ1,{'all'},3);
%view([0,90])

OBJT2 = transformOBJ(OBJ1,-pi/2,[0,0,1]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(OutputO.cleanX,-pi/2,[0,0,1]);
    OutputO.cleanX = Quaternion3(-pi/2,[0,0,1],OutputO.cleanX);
end
%h = PlotEmbryoCS7GP(Output,OBJT2,{'all'},2);

OBJT3 = transformOBJ(OBJT2,1.2*pi,[0,1,0]);
if length(varargin)==3
%OutputO.cleanX = transformOBJ(OutputO.cleanX,1.2*pi,[0,1,0]);
    OutputO.cleanX = Quaternion3(1.2*pi,[0,1,0],OutputO.cleanX);

end
%h = PlotEmbryoCS7GP(Output,OBJT3,{'all'},3);
%view([-7.5153,90])

aa = linspace(0,2*pi,100);
i = 15%:100
OBJout = transformOBJ(OBJT3,aa(i),[1,0,.1]/norm([1,0,.1]));
if length(varargin)==3
%OutputO.cleanX = transformOBJ(OutputO.cleanX,aa(i),[1,0,.1]/norm([1,0,.1]));
    OutputO.cleanX = Quaternion3(aa(i),[1,0,.1]/norm([1,0,.1]),OutputO.cleanX);

end
%h2 = PlotEmbryoCS7GP(Output,OBJT4,{'all'},4);
%view([-7.5153,90])
%camlight('left')
%pause
%clf
%end
a = -11.0469;
b = 60.5571;
else
%h2 = PlotEmbryoCS7GP(Output,OBJ1,{'EmDisc','Stalk'},4);

%view([-3.1291,22.6770])
%view([-1.9291,36.0849])
a = -1.0233;
b = 33.8659; %]) %+flip x
OBJout = OBJ1;
%Whole disc viiew 
%a = -33.3166; 
%b = 36.2020;

%Whole disc viiew [-33.3166,36.2020]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end