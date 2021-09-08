function [OBJout,a,b,OutputO] = transformCS5VE(varargin)

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

%View for CS5 VE only
OBJ1c_2 = transformOBJ(OBJ1,'flipZ')
OBJout = transformOBJ(OBJ1c_2,'flipX')


% h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c_2,{'VE'},1);
% view([138.5913,5.2765])
% %camlight('left')
% view([c1,d1])
% ax = gca(h2)
% ax.CLim = [mi1 mi2];
% view([-183.7337,-15.1533])
% camlight('left')
% view([-179.8734,-61.7588])
a = -179.8734;
b = -61.7588;

else

%View for VE and EmDisc CS5
OBJ1c_3 = transformOBJ(OBJ1,'flipX')
% h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c_3,{'VE','EmDisc'},4);
% view([138.5913,5.2765])
% %camlight('left')
% view([c1,d1])
% ax = gca(h2)
% ax.CLim = [mi1 mi2];
% view([-183.7337,-15.1533])
% view([-179.8734,-61.7588])
% view([8.6708,30.1752])
% camlight('right')
%view([
a = 4.5723;
b = 77.6444;

end

% 
% 
% 
% 
% 
% aa = linspace(0,-2*pi,100);
% 
% i = 11;% %1:100
% OBJT2 = transformOBJ(OBJ1,aa(i),[1,0,0]);
% if length(varargin)==3
%     OutputO.cleanX = Quaternion3(aa(i),[1,0,0],OutputO.cleanX);
% end
% %h = PlotEmbryoCS5GP(Output,OBJT2,{'all'},2);%
% %view([183.0698,90])
% %pause
% %clf
% %end
% %camlight('right')
% %camlight('left')
% %i = 11;
% %view([162.4898,83.3904])
% %view([177.9950,80.9869])
% 
% OBJout = transformOBJ(OBJT2,'flipX');
% if length(varargin)==3
% OutputO.cleanX(:,1) = -OutputO.cleanX(:,1);
% end
% %h = PlotEmbryoCS5GP(Output,OBJT3,{'all'},3);
% %view([177.9950,80.9869])
% %camlight('left')
% %view([
% a = 164.8625;
% b = 74.9781;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% else
% aa = linspace(0,-2*pi,100);
% 
% i = 11; %1:100
% OBJT2 = transformOBJ(OBJ1,aa(i),[1,0,0]);
% if length(varargin)==3
%     OutputO.cleanX = Quaternion3(aa(i),[1,0,0],OutputO.cleanX);
% end
% %h = PlotEmbryoCS5GP(Output,OBJT2,{'EmDisc','VE'},2);
% %view([183.0698,90])
% %pause
% %clf
% %end
% %camlight('right')
% %camlight('left')
% %i = 11;
% %view([162.4898,83.3904])
% %view([177.9950,80.9869])
% 
% OBJout = transformOBJ(OBJT2,'flipZ');
% if length(varargin)==3
% OutputO.cleanX(:,3) = -OutputO.cleanX(:,3);
% end
% %h = PlotEmbryoCS5GP(Output,OBJT3,{'EmDisc','VE'},3);%
% %view([177.9950,80.9869])
% %camlight('left')
% %view([7.1060,81.5877])
% 
% %[Output] = MarmosetGP_CS5(D,Output,'MIXL1');
% %[Output] = MarmosetGPInfer_CS5(Output,OBJ1);
% %h = PlotEmbryoCS5GP(Output,OBJT3,{'EmDisc','VE'},4);
% %view([7.1060,81.5877])
% %camlight('left')
% %view([
% a = -5.9798;
% b = 89.3991;
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%
% % aa = linspace(0,-2*pi,100);
% % for i = 1:100
% % OBJT3 = transformOBJ(OBJT2,aa(i),[0,1,0]);
% % h = PlotEmbryoCS5GP(Output,OBJT3,{'all'},3);
% % view([177.9950,80.9869])
% % camlight('left')
% % pause
% % clf
% % end
% % 
% % 
% % aa = linspace(0,-2*pi,100);
% % for i = 1:100
% % OBJT2 = transformOBJCS5(OBJ1,aa(i),[0,1,0]);
% % h = PlotEmbryoCS5GP(Output,OBJT2,{'all'},2);
% % view([183.0698,90])
% % pause
% % clf
% % end
% % camlight('right')
% % camlight('left')
% % 
% % 
% % 
% % %view([1.4946,-85.7939])
% % 
% % % 183.0698,90
% %  
% % %view([183.5224,90])
% % 
% % view([-0.0107,-0.0454,5.1247])
% end