addpath(genpath('./Functions'))

%Load CS7
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,Shots);
%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS6');
[OutputCS6] = loadCS6Scaffold(D,Locations,Shots);
%Load CS6
[OBJ2_2,section] = LoadCS6_2('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS62');
[OutputCS62] = loadCS6Scaffold2(D,Locations,Shots);
%Load CS7
[OBJ3,section] = LoadCS7('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS7');
[OutputCS7] = loadCS7Scaffold(D,Locations,Shots);
OutputCS7.scalefactor = 600;

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ3b,a3,b3,OutputCS7b] = transformCS7(OBJ3,OutputCS7,'all');
[OBJ2_2b,a2_2,b2_2,OutputCS62b] = transformCS62(OBJ2_2,OutputCS62,'all');

[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');
[OBJ3c,c3,d3,OutputCS7c] = transformCS7(OBJ3,OutputCS7,'notall');
[OBJ2_2c,c2_2,d2_2,OutputCS62c] = transformCS62(OBJ2_2,OutputCS62,'notall');

gene = {'HHEX'}


%First do standard run scale by 1200
for i = 1:length(gene)
    
    try

[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,gene{i});
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);

HYPCS{i,1} = OutputCS5.HYP1;
HYPCS{i,2} = OutputCS5.HYP2;

L1{i,1} = OutputCS5.L1;
L1{i,2} = OutputCS5.L2;

[OutputCS6] = MarmosetGP_CS6_v3(D,OutputCS6,gene{i});
[OutputCS6] = MarmosetGPInfer_CS6_v3(OutputCS6,OBJ2);

HYPCS{i,3} = OutputCS6.HYP1;
HYPCS{i,4} = OutputCS6.HYP2;

L1{i,3} = OutputCS6.L1;
L1{i,4} = OutputCS6.L2;

[OutputCS7] = MarmosetGP_CS7_v3(D,OutputCS7,gene{i});
[OutputCS7] = MarmosetGPInfer_CS7_v3(OutputCS7,OBJ3);

L1{i,5} = OutputCS7.L1;
L1{i,6} = OutputCS7.L2;

HYPCS{i,5} = OutputCS7.HYP1;
HYPCS{i,6} = OutputCS7.HYP2;

[OutputCS62] = MarmosetGP_CS62(D,OutputCS62,gene{i});
[OutputCS62] = MarmosetGPInfer_CS6(OutputCS62,OBJ2_2);

HYPCS{i,7} = OutputCS62.HYP1;
HYPCS{i,8} = OutputCS62.HYP2;

L1{i,7} = OutputCS62.L1;
L1{i,8} = OutputCS62.L2;

m1 = OutputCS5.cLim;
m2 = OutputCS6.cLim;
m3 = OutputCS7.cLim;

ma1s = max(max([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));
mi1s = min(min([OutputCS5.cLim2;OutputCS5.cLim4;OutputCS6.cLim2;OutputCS6.cLim2b;OutputCS6.cLim4;OutputCS7.cLim2]));

mi1 = min([m1,m2,m3]);
ma1 = max([m1,m2,m3]);

m2 = max([ma1,0.07]);
m1 = max([mi1,0]);
mi2 = max([ma1s,0.07]);
mi1 = max([mi1s,0]);

h1 = PlotEmbryoCS5GP(OutputCS5,OBJ1b,{'all'},2,[5,8,1,2,9,10]);
view(a1,b1)
camlight('left')
ax = gca(h1)
ax.CLim = [m1 m2];

h3 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2b,{'all'},1);
view([-34.0286,-68.0304])
camlight('headlight')
view(a2,b2)
ax = gca(h3)
ax.XLim = ax.XLim*0.9;
ax.YLim = ax.YLim*0.9;
ax.CLim = [m1 m2];

h3 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2b,{'all'},2,[5,8,3,4,11,12]);
view([269.5795,2.8730])
camlight('left')
view(a2,b2)
ax = gca(h3)
ax.XLim = ax.XLim*0.9;
ax.YLim = ax.YLim*0.9;
ax.CLim = [m1 m2];

h5 = PlotEmbryoCS7GP_v3(OutputCS7,OBJ3b,{'all'},2,[5,8,5,6,13,14]);
view(a3,b3)
camlight('left')
ax = gca(h5)
ax.XLim = ax.XLim*0.9;
ax.YLim = ax.YLim*0.9;
ax.CLim = [m1 m2];
cb=colorbar;
cb.Position = cb.Position + [0.03 0.05 -0.002 -0.1];
cb.Color = 'k';
cb.FontSize = 12;
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(end)])
end

h3b = PlotEmbryoCS6GP_2(OutputCS62,OBJ2_2b,{'all'},2,[5,8,7,8,15,16]);
view(a2_2,b2_2)
camlight('left')
ax = gca(h3b)
ax.CLim = [m1 m2];

OBJ1c_2 = transformOBJ(OBJ1c,'flipZ')
OBJ1c_2 = transformOBJ(OBJ1c_2,'flipX')
OBJ1c_3 = transformOBJ(OBJ1c_2,'flipX')
h4 = PlotEmbryoCS5GP(OutputCS5,OBJ1c_3,{'VE','EmDisc'},2,[5,8,17,18,25,26]);
ax = gca(h4)
ax.CLim = [mi1 mi2];
view([8.6708,30.1752])
camlight('right')
view([4.5723,77.6444])
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;

h4 = PlotEmbryoCS6GP(OutputCS6,OBJ2c,{'VE','Stalk','EmDisc'},2,[5,8,19,20,27,28]);
ax = gca(h4)
ax.CLim = [mi1 mi2];
view([ -86.0444,-90])
camlight('left')
ax.XLim = ax.XLim*1.4;
ax.YLim = ax.YLim*1.4;

h6 = PlotEmbryoCS7GP_v3(OutputCS7,OBJ3c,{'EmDisc','Stalk'},2,[5,8,21,22,29,30]);
view([77.0455, -27.1350])
camlight('left')
view([c3,d3])
ax = gca(h6)
ax.CLim = [mi1 mi2];
cb=colorbar;
cb.Position = cb.Position + [0.03 0.05 -0.002 -0.1];
cb.Color = 'k';
cb.FontSize = 12;
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(end)])
end

%CS6 EmDisc, Stalk, VE
OBJ2_2c_2 = transformOBJ(OBJ2_2c,'flipY')
OBJ2_2c_2 = transformOBJ(OBJ2_2c_2,0.6*pi,[1,0,0])
h4b = PlotEmbryoCS6GP_2_v3(OutputCS62,OBJ2_2c_2,{'EmDisc','VE','Stalk'},2,[5,8,23,24,31,32]);
ax = gca(h4b)
ax.CLim = [mi1 mi2];
view([180.045,180.0459])
camlight('right')
view([124.1610,79.7851])

load('../Data/SpatialData/CS5_VE.mat')
[Output5] = MarmosetGPInfer_CS5(OutputCS5,Line,'Line');
ha1 = plotAP(Output5, 2, [5,8,33,34], 'VE');
axis square
pbaspect([1.5 1 1])


load('../Data/SpatialData/CS6_VE.mat')
[Output6] = MarmosetGPInfer_CS6_v3(OutputCS6,Line,'Line');
ha2 = plotAP(Output6, 2, [5,8,35,36], 'VE');
%axis square
pbaspect([1.5 1 1])

linkaxes([ha1 ha2],'xy')
ha1.YLim = ha1.YLim*1.1;

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
print('-dpng',['./Plots/' gene{i} '_APVE.png'],'-r1000')

clf

load('../Data/SpatialData/CS5_VE.mat')
[Output5] = MarmosetGPInfer_CS5(OutputCS5,Line,'Line');
ha1 = plotAP(Output5, 2, [5,8,33,34], 'VE');
axis square
pbaspect([1.5 1 1])


load('../Data/SpatialData/CS6_VE.mat')
[Output6] = MarmosetGPInfer_CS6_v3(OutputCS6,Line,'Line');
ha2 = plotAP(Output6, 2, [5,8,35,36], 'VE');
%axis square
pbaspect([1.5 1 1])

linkaxes([ha1 ha2],'xy')
ha1.YLim = ha1.YLim*1.1;

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
print('-dpdf',['./Plots/' gene{i} '_APVVELP.pdf'],'-r1000')

close all

    end
end