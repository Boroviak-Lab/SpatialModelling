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

%gene = {'SOX2','POU5F1','PRDM14','T','LEF1','SFRP1','SFRP2','SFRP5','NANOG','WNT3','WNT8A','FBXO2','TDGF1','MIXL1','GATA6','SOX17','CDH2','CER1','LEFTY2','LHX1','OTX2','NODAL','SDC4','PDGFRA','PDGFA','KDR','VEGFA'};
gene = unique({'SOX2','NANOG','POU5F1','MIXL1','T','EOMES','GATA6','LEF1','PDGFRA','ID2','NODAL','CER1','LEFTY2','LHX1','OTX2','HHEX','NANOS3','TFAP2A','TFAP2C','VTCN1','TTR','GC','GATA2','GATA3','HAND2','HGF','ID1','ID3','WNT5A','WT5B','WNT6','TCF4','DKK1','ID1','ID2','NOG','BAMBI','BMP2','BMP6','SOX17','APOA1','GC','ARL13B','IHH','PCH2','SMO','GPC4','GLI1','SOX17','KLF4','DAZL','MAEL','PRAME','NODAL','TDGF1','BMP4','ID2','SOX2','POU5F1','PRDM14','T','LEF1','SFRP1','SFRP2','SFRP5','NANOG','WNT3','WNT8A','FBXO2','TDGF1','MIXL1','GATA6','SOX17','CDH2','CER1','LEFTY2','LHX1','OTX2','NODAL','SDC4','PDGFRA','PDGFA','KDR','VEGFA'});

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ3b,a3,b3,OutputCS7b] = transformCS7(OBJ3,OutputCS7,'all');
[OBJ2_2b,a2_2,b2_2,OutputCS62b] = transformCS62(OBJ2_2,OutputCS62,'all');

[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');
[OBJ3c,c3,d3,OutputCS7c] = transformCS7(OBJ3,OutputCS7,'notall');
[OBJ2_2c,c2_2,d2_2,OutputCS62c] = transformCS62(OBJ2_2,OutputCS62,'notall');

%gene = {'HHEX'}


%1 2 3 4        5 6 7 8     9 10 11 12   13 14 15 16
%17 18 19 20   21 22 23 24  25 26 27 28  29 30 31 32
%33 34 35 36 . 37 38 39 40 . 41 42 43 44 . 45 46 47 48
%49 50 51 52   53 54 55 56 . 57 58 59 60 . 61 62 63 64



%First do standard run scale by 1200
for i = 1:length(gene)
    
    try

[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,gene{i});
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);

HYPCS{i,1} = OutputCS5.HYP1;
HYPCS{i,2} = OutputCS5.HYP2;

L1{i,1} = OutputCS5.L1;
L1{i,2} = OutputCS5.L2;

[OutputCS6] = MarmosetGP_CS6(D,OutputCS6,gene{i});
[OutputCS6] = MarmosetGPInfer_CS6_v3(OutputCS6,OBJ2);

HYPCS{i,3} = OutputCS6.HYP1;
HYPCS{i,4} = OutputCS6.HYP2;

L1{i,3} = OutputCS6.L1;
L1{i,4} = OutputCS6.L2;

[OutputCS7] = MarmosetGP_CS7(D,OutputCS7,gene{i});
[OutputCS7] = MarmosetGPInfer_CS7(OutputCS7,OBJ3);

L1{i,5} = OutputCS7.L1;
L1{i,6} = OutputCS7.L2;

HYPCS{i,5} = OutputCS7.HYP1;
HYPCS{i,6} = OutputCS7.HYP2;

[OutputCS62] = MarmosetGP_CS62(D,OutputCS62,gene{i});
[OutputCS62] = MarmosetGPInfer_CS6_v3(OutputCS62,OBJ2_2);

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
%ax.CLim = [max(OutputCS5.cLim(1),0) max(OutputCS5.cLim(1),0.07)];
ax.CLim = [m1 m2];


h3 = PlotEmbryoCS6GP(OutputCS6,OBJ2b,{'all'},2,[5,8,3,4,11,12]);
view(a2,b2)
camlight('left')
ax = gca(h3)
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];




h5 = PlotEmbryoCS7GP(OutputCS7,OBJ3b,{'all'},2,[5,8,5,6,13,14]);
view(a3,b3)
camlight('left')
ax = gca(h5)
%ax.CLim = [max(OutputCS7.cLim(1),0) max(OutputCS7.cLim(1),0.07)];
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
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];


%View for CS5 VE only
%OBJ1c_2 = transformOBJ(OBJ1c,'flipZ')
%OBJ1c_2 = transformOBJ(OBJ1c_2,'flipX')
%h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c_2,{'VE'},1);
%view([138.5913,5.2765])
%%camlight('left')
%view([c1,d1])
%ax = gca(h2)
%ax.CLim = [mi1 mi2];
%view([-183.7337,-15.1533])
%camlight('left')
%view([-179.8734,-61.7588])

%%View for VE and EmDisc CS5
%OBJ1c_3 = transformOBJ(OBJ1c_2,'flipX')
%h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c_3,{'VE','EmDisc'},4);
%view([138.5913,5.2765])
%%camlight('left')
%view([c1,d1])
%ax = gca(h2)
%ax.CLim = [mi1 mi2];
%view([-183.7337,-15.1533])
%view([-179.8734,-61.7588])
%view([8.6708,30.1752])
%camlight('right')
%view([4.5723,77.6444])%




h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c,{'EmDisc','VE'},2,[5,8,17,18,25,26]);
view([138.5913,5.2765])
camlight('left')
view([c1,d1])
ax = gca(h2)
ax.CLim = [mi1 mi2];
ax.XLim = ax.XLim*1.2;
ax.YLim = ax.YLim*1.2;


h4 = PlotEmbryoCS6GP(OutputCS6,OBJ2c,{'EmDisc','VE','Stalk'},2,[5,8,19,20,27,28]);
view([8.5335,38.4802])
camlight('left')
view([c2,d2])
ax = gca(h4)
ax.CLim = [mi1 mi2];
ax.XLim = ax.XLim*1.4;
ax.YLim = ax.YLim*1.4;
axis on

%CS6 VE,Stalk, EmD 
%h4 = PlotEmbryoCS6GP(OutputCS6,OBJ2c,{'VE','Stalk','EmDisc'},1);
%view([8.5335,38.4802])
%%camlight('left')
%view([c2,d2])
%ax = gca(h4)
%ax.CLim = [mi1 mi2];
%view([ -86.0444,-90])
%camlight('left')

%h4 = PlotEmbryoCS6GP(OutputCS6,OBJ2c,{'VE'},1);
%view([8.5335,38.4802])
%%camlight('left')
%view([c2,d2])
%ax = gca(h4)
%ax.CLim = [mi1 mi2];
%view([ -86.0444,-90])
%camlight('left')
%view([-24.6146,83.6531])



h6 = PlotEmbryoCS7GP(OutputCS7,OBJ3c,{'EmDisc','Stalk'},2,[5,8,21,22,29,30]);
%view([84.9000,28.6748])
%view([92.3512,20.2910])
 %view([111.3449,67.4202])
view([77.0455, -27.1350])
camlight('left')
view([c3,d3])
%camlight('left')
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




h4b = PlotEmbryoCS6GP_2(OutputCS62,OBJ2_2c,{'EmDisc','VE','Stalk'},2,[5,8,23,24,31,32]);
view([c2_2,d2_2])
camlight('left')
ax = gca(h4b)
ax.CLim = [mi1 mi2];

%CS6 EmDisc, Stalk, VE
%OBJ2_2c_2 = transformOBJ(OBJ2_2c,'flipY')
%OBJ2_2c_2 = transformOBJ(OBJ2_2c_2,0.6*pi,[1,0,0])
%h4b = PlotEmbryoCS6GP_2(OutputCS62,OBJ2_2c_2,{'EmDisc','VE','Stalk'},1);
%view([c2_2,d2_2])
%%camlight('left')
%%ax = gca(h4b)
%ax.CLim = [mi1 mi2];
%view([180.045,180.0459])
%%view([248.6674,-17.5568])
%camlight('right')%
%view([124.1610,79.7851])

%CS6_2 Ve only
%OBJ2_2c_2 = transformOBJ(OBJ2_2c,'flipY')
%OBJ2_2c_2 = transformOBJ(OBJ2_2c_2,0.6*pi,[1,0,0])
%%h4b = PlotEmbryoCS6GP_2(OutputCS62,OBJ2_2c_2,{'VE'},1);
%view([c2_2,d2_2])
%%camlight('left')
%ax = gca(h4b)
%ax.CLim = [mi1 mi2];
%view([180.045,180.0459])
%view([248.6674,-17.5568])
%view([124.1610,79.7851])
%view([45.4144,-90])
%camlight('right')

load('./Data/CS5_EmDisc.mat')
[Output5] = MarmosetGPInfer_CS5(OutputCS5,Line,'Line');
ha1 = plotAP(Output5, 2, [5,8,33,34]);
%[5,8,33,34]
%[5,8,35,36]
%[5,8,37,38]
%[5,8,39,40]

load('./Data/CS6_EmDisc.mat')
[Output6] = MarmosetGPInfer_CS6(OutputCS6,Line,'Line');
ha2 = plotAP(Output6, 2, [5,8,35,36]);

load('./Data/CS7_EmDisc.mat')
[Output7] = MarmosetGPInfer_CS7(OutputCS7,Line,'Line');
ha3 = plotAP(Output7, 2, [5,8,37,38]);
linkaxes([ha1 ha2 ha3],'xy')
ha1.YLim = ha1.YLim*1.1;
%5 x 8
%1 2   .  3 4     5 6       7 8
%9 10 .  11 12 .  13 14 . 15 16
%17 18 . 19 20 . 21 22 . 23 24
%25 26 . 27 28 . 29 30 . 31 32
%33 34 . 35 36 . 37 38 . 39 40


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
print('-dpng',['~/Desktop/Test1.png'],'-r1000')

%print(h1,['Plots/WholeEmb_' gene{i} '_CS5_scaled1200.pdf'],'-dpdf','-r1000');
%print(h3,['Plots/WholeEmb_' gene{i} '_CS6_scaled1200.pdf'],'-dpdf','-r1000');
print(h5,['Plots/WholeEmb_' gene{i} '_CS5-7_scaled600.pdf'],'-dpdf','-r1000');

%print(h2,['Plots/AP_' gene{i} '_CS5_scaled1200.pdf'],'-dpdf','-r1000');
%print(h4,['Plots/AP_' gene{i} '_CS6_scaled1200.pdf'],'-dpdf','-r1000');
print(h6,['Plots/AP_' gene{i} '_CS5-7_scaled600.pdf'],'-dpdf','-r1000');

close all

    end
end

save('FortyGene_Stalkmerged_scale600_L1.mat','L1')
save('FortyGene_Stalkmerged_scale600_Hyp.mat','HYPCS')


%Now standard run but scale by 800
OutputCS7.scalefactor = 600;
for i = 1:length(gene)
    
    try

[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,gene{i});
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);

HYPCS{i,1} = OutputCS5.HYP1;
HYPCS{i,2} = OutputCS5.HYP2;

L1{i,1} = OutputCS5.L1;
L1{i,2} = OutputCS5.L2;

[OutputCS6] = MarmosetGP_CS6(D,OutputCS6,gene{i});
[OutputCS6] = MarmosetGPInfer_CS6(OutputCS6,OBJ2);

HYPCS{i,3} = OutputCS6.HYP1;
HYPCS{i,4} = OutputCS6.HYP2;

L1{i,3} = OutputCS6.L1;
L1{i,4} = OutputCS6.L2;

[OutputCS7] = MarmosetGP_CS7(D,OutputCS7,gene{i});
[OutputCS7] = MarmosetGPInfer_CS7(OutputCS7,OBJ3);

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


h1 = PlotEmbryoCS5GP(OutputCS5,OBJ1b,{'all'},1,[2,2,1]);
view(a1,b1)
camlight('left')
ax = gca(h1)
%ax.CLim = [max(OutputCS5.cLim(1),0) max(OutputCS5.cLim(1),0.07)];
ax.CLim = [m1 m2];

h3 = PlotEmbryoCS6GP(OutputCS6,OBJ2b,{'all'},1,[2,2,2]);
view(a2,b2)
camlight('left')
ax = gca(h3)
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];


h5 = PlotEmbryoCS7GP(OutputCS7,OBJ3b,{'all'},1,[2,2,3]);
view(a3,b3)
camlight('left')
ax = gca(h5)
%ax.CLim = [max(OutputCS7.cLim(1),0) max(OutputCS7.cLim(1),0.07)];
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


h3b = PlotEmbryoCS6GP_2(OutputCS62,OBJ2_2b,{'all'},1,[2,2,4]);
view(a2_2,b2_2)
camlight('left')
ax = gca(h3b)
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];


h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c,{'EmDisc','VE'},2,[2,2,1]);
view([138.5913,5.2765])
camlight('left')
view([c1,d1])
ax = gca(h2)
ax.CLim = [mi1 mi2];

h4 = PlotEmbryoCS6GP(OutputCS6,OBJ2c,{'EmDisc','VE','Stalk'},2,[2,2,2]);
view([8.5335,38.4802])
camlight('left')
view([c2,d2])
ax = gca(h4)
ax.CLim = [mi1 mi2];


h6 = PlotEmbryoCS7GP(OutputCS7,OBJ3c,{'EmDisc','Stalk'},2,[2,2,3]);
%view([84.9000,28.6748])
%view([92.3512,20.2910])
 %view([111.3449,67.4202])
view([77.0455, -27.1350])
camlight('left')
view([c3,d3])
%camlight('left')
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

h4b = PlotEmbryoCS6GP_2(OutputCS62,OBJ2_2c,{'EmDisc','VE','Stalk'},2,[2,2,4]);
view([c2_2,d2_2])
camlight('left')
ax = gca(h4)
ax.CLim = [mi1 mi2];


%print(h1,['Plots/WholeEmb_' gene{i} '_CS5_scaled800.pdf'],'-dpdf','-r1000');
%print(h3,['Plots/WholeEmb_' gene{i} '_CS6_scaled800.pdf'],'-dpdf','-r1000');
print(h5,['Plots/WholeEmb_' gene{i} '_CS5-7_scaled800.pdf'],'-dpdf','-r1000');

%print(h2,['Plots/AP_' gene{i} '_CS5_scaled800.pdf'],'-dpdf','-r1000');
%print(h4,['Plots/AP_' gene{i} '_CS6_scaled800.pdf'],'-dpdf','-r1000');
print(h6,['Plots/AP_' gene{i} '_CS5-7_scaled800.pdf'],'-dpdf','-r1000');

close all

    end
end

save('FortyGene_Stalkmerged_scale800_L1.mat','L1')
save('FortyGene_Stalkmerged_scale800_Hyp.mat','HYPCS')


%Now seperate out the stalk. Scale by 800
OutputCS7.scalefactor = 600;
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

[OutputCS62] = MarmosetGP_CS62_v3(D,OutputCS62,gene{i});
[OutputCS62] = MarmosetGPInfer_CS6_v3(OutputCS62,OBJ2_2);

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


h1 = PlotEmbryoCS5GP(OutputCS5,OBJ1b,{'all'},1,[2,2,1]);
view(a1,b1)
camlight('left')
ax = gca(h1)
%ax.CLim = [max(OutputCS5.cLim(1),0) max(OutputCS5.cLim(1),0.07)];
ax.CLim = [m1 m2];

h3 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2b,{'all'},1,[2,2,2]);
view(a2,b2)
camlight('left')
ax = gca(h3)
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];


h5 = PlotEmbryoCS7GP_v3(OutputCS7,OBJ3b,{'all'},1,[2,2,3]);
view(a3,b3)
camlight('left')
ax = gca(h5)
%ax.CLim = [max(OutputCS7.cLim(1),0) max(OutputCS7.cLim(1),0.07)];
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


h3b = PlotEmbryoCS6GP_2_v3(OutputCS62,OBJ2_2b,{'all'},1,[2,2,4]);
view(a2_2,b2_2)
camlight('left')
ax = gca(h3b)
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];


h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c,{'EmDisc','VE'},2,[2,2,1]);
view([138.5913,5.2765])
camlight('left')
view([c1,d1])
ax = gca(h2)
ax.CLim = [mi1 mi2];

h4 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2c,{'EmDisc','VE','Stalk'},2,[2,2,2]);
view([8.5335,38.4802])
camlight('left')
view([c2,d2])
ax = gca(h4)
ax.CLim = [mi1 mi2];


h6 = PlotEmbryoCS7GP_v3(OutputCS7,OBJ3c,{'EmDisc','Stalk'},2,[2,2,3]);
%view([84.9000,28.6748])
%view([92.3512,20.2910])
 %view([111.3449,67.4202])
view([77.0455, -27.1350])
camlight('left')
view([c3,d3])
%camlight('left')
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

h4b = PlotEmbryoCS6GP_2_v3(OutputCS62,OBJ2_2c,{'EmDisc','VE','Stalk'},2,[2,2,4]);
view([c2_2,d2_2])
camlight('left')
ax = gca(h4b)
ax.CLim = [mi1 mi2];


%print(h1,['Plots/WholeEmb_' gene{i} '_CS5_scaled800.pdf'],'-dpdf','-r1000');
%print(h3,['Plots/WholeEmb_' gene{i} '_CS6_scaled800.pdf'],'-dpdf','-r1000');
print(h5,['Plots/WholeEmb_' gene{i} '_CS5-7_sep600.pdf'],'-dpdf','-r1000');

%print(h2,['Plots/AP_' gene{i} '_CS5_scaled800.pdf'],'-dpdf','-r1000');
%print(h4,['Plots/AP_' gene{i} '_CS6_scaled800.pdf'],'-dpdf','-r1000');
print(h6,['Plots/AP_' gene{i} '_CS5-7_sep600.pdf'],'-dpdf','-r1000');

close all

    end
end

save('FortyGene_Stalksep_scale600_L1.mat','L1')
save('FortyGene_Stalksep_scale600_Hyp.mat','HYPCS')


%Now scale by 1200
OutputCS7.scalefactor = 1200;
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

[OutputCS62] = MarmosetGP_CS62_v3(D,OutputCS62,gene{i});
[OutputCS62] = MarmosetGPInfer_CS6_v3(OutputCS62,OBJ2_2);

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


h1 = PlotEmbryoCS5GP(OutputCS5,OBJ1b,{'all'},1,[2,2,1]);
view(a1,b1)
camlight('left')
ax = gca(h1)
%ax.CLim = [max(OutputCS5.cLim(1),0) max(OutputCS5.cLim(1),0.07)];
ax.CLim = [m1 m2];

h3 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2b,{'all'},1,[2,2,2]);
view(a2,b2)
camlight('left')
ax = gca(h3)
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];


h5 = PlotEmbryoCS7GP_v3(OutputCS7,OBJ3b,{'all'},1,[2,2,3]);
view(a3,b3)
camlight('left')
ax = gca(h5)
%ax.CLim = [max(OutputCS7.cLim(1),0) max(OutputCS7.cLim(1),0.07)];
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


h3b = PlotEmbryoCS6GP_2_v3(OutputCS62,OBJ2_2b,{'all'},1,[2,2,4]);
view(a2_2,b2_2)
camlight('left')
ax = gca(h3b)
%ax.CLim = [max(OutputCS6.cLim(1),0) max(OutputCS6.cLim(1),0.07)];
ax.CLim = [m1 m2];


h2 = PlotEmbryoCS5GP(OutputCS5,OBJ1c,{'EmDisc','VE'},2,[2,2,1]);
view([138.5913,5.2765])
camlight('left')
view([c1,d1])
ax = gca(h2)
ax.CLim = [mi1 mi2];

h4 = PlotEmbryoCS6GP_v3(OutputCS6,OBJ2c,{'EmDisc','VE','Stalk'},2,[2,2,2]);
view([8.5335,38.4802])
camlight('left')
view([c2,d2])
ax = gca(h4)
ax.CLim = [mi1 mi2];


h6 = PlotEmbryoCS7GP_v3(OutputCS7,OBJ3c,{'EmDisc','Stalk'},2,[2,2,3]);
%view([84.9000,28.6748])
%view([92.3512,20.2910])
 %view([111.3449,67.4202])
view([77.0455, -27.1350])
camlight('left')
view([c3,d3])
%camlight('left')
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

h4b = PlotEmbryoCS6GP_2_v3(OutputCS62,OBJ2_2c,{'EmDisc','VE','Stalk'},2,[2,2,4]);
view([c2_2,d2_2])
camlight('left')
ax = gca(h4)
ax.CLim = [mi1 mi2];


%print(h1,['Plots/WholeEmb_' gene{i} '_CS5_scaled800.pdf'],'-dpdf','-r1000');
%print(h3,['Plots/WholeEmb_' gene{i} '_CS6_scaled800.pdf'],'-dpdf','-r1000');
print(h5,['Plots/WholeEmb_' gene{i} '_CS5-7_sep1200.pdf'],'-dpdf','-r1000');

%print(h2,['Plots/AP_' gene{i} '_CS5_scaled800.pdf'],'-dpdf','-r1000');
%print(h4,['Plots/AP_' gene{i} '_CS6_scaled800.pdf'],'-dpdf','-r1000');
print(h6,['Plots/AP_' gene{i} '_CS5-7_sep1200.pdf'],'-dpdf','-r1000');

close all

    end
end

save('FortyGene_Stalksep_sep1200_L1.mat','L1')
save('FortyGene_Stalksep_sep1200_Hyp.mat','HYPCS')



[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ3b,a3,b3,OutputCS7b] = transformCS7(OBJ3,OutputCS7,'all');
[OBJ2_2b,a2_2,b2_2,OutputCS62b] = transformCS62(OBJ2_2,OutputCS62,'all');

[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');
[OBJ3c,c3,d3,OutputCS7c] = transformCS7(OBJ3,OutputCS7,'notall');
[OBJ2_2c,c2_2,d2_2,OutputCS62c] = transformCS62(OBJ2_2,OutputCS62,'notall');


%Now do scatterplots
OutputCS7.scalefactor = 1200;
for i = 1:length(gene)
    
    try

[OutputCS5] = MarmosetGP_CS5(D,OutputCS5,gene{i});
[OutputCS5] = MarmosetGPInfer_CS5(OutputCS5,OBJ1);

[OutputCS6] = MarmosetGP_CS6_v3(D,OutputCS6,gene{i});
[OutputCS6] = MarmosetGPInfer_CS6_v3(OutputCS6,OBJ2);

[OutputCS7] = MarmosetGP_CS7_v3(D,OutputCS7,gene{i});
[OutputCS7] = MarmosetGPInfer_CS7_v3(OutputCS7,OBJ3);

[OutputCS62] = MarmosetGP_CS62_v3(D,OutputCS62,gene{i});
[OutputCS62] = MarmosetGPInfer_CS6_v3(OutputCS62,OBJ2_2);

m1 = min([OutputCS5.Ytrain,OutputCS6.Ytrain,OutputCS7.Ytrain]);
m2 = max([OutputCS5.Ytrain,OutputCS6.Ytrain,OutputCS7.Ytrain]);

hh = figure(1)
h1 = subplot(2,2,1);
scatter3(OutputCS5b.cleanX(:,1),OutputCS5b.cleanX(:,2),OutputCS5b.cleanX(:,3),125,OutputCS5.Ytrain','filled')
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(a1,b1)



%        indT = find(strcmp(Output.cleanAnotaton,'Am_CS5')==1 | strcmp(Output.cleanAnotaton,'Am_CS5_PGC')==1);
%        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,'ko','filled')
%        scatter3(XYZ(indT,1),XYZ(indT,2), XYZ(indT,3),125,Y(indT),'filled')


h2 = subplot(2,2,2);
scatter3(OutputCS6b.cleanX(:,1),OutputCS6b.cleanX(:,2),OutputCS6b.cleanX(:,3),125,OutputCS6.Ytrain,'filled')
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(a2,b2)

subplot(2,2,3);
scatter3(OutputCS7b.cleanX(:,1),OutputCS7b.cleanX(:,2),OutputCS7b.cleanX(:,3),125,OutputCS7.Ytrain,'filled')
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(a3,b3)

subplot(2,2,4);
scatter3(OutputCS62b.cleanX(:,1),OutputCS62b.cleanX(:,2),OutputCS62b.cleanX(:,3),125,OutputCS62.Ytrain,'filled')
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(a2_2,b2_2)


print(hh,['Plots/WholeEmbScatter_' gene{i} '_CS5-7.pdf'],'-dpdf','-r1000');

hh = figure(2);

indT1 =  find(strcmp(OutputCS5.cleanAnotaton,'EmDisc_CS5')==1  | strcmp(OutputCS5.cleanAnotaton,'VE_CS5')==1);
indT2 = find(strcmp(OutputCS6.cleanAnotaton,'EmDisc_CS6')==1 | strcmp(OutputCS6.cleanAnotaton,'VE_CS6')==1 | strcmp(OutputCS6.cleanAnotaton,'EmDisc_Gast_CS6')==1 | strcmp(OutputCS6.cleanAnotaton,'EmDisc_Stalk_CS6')==1 | strcmp(OutputCS6.cleanAnotaton,'EmDisc_gast_CS6')==1 | strcmp(OutputCS6.cleanAnotaton,'EmDisc_stalk_CS6')==1 | strcmp(OutputCS6.cleanAnotaton,'VE_CS6')==1);
indT3 = find(strcmp(OutputCS7.cleanAnotaton,'EmDisc_CS7')==1 | strcmp(OutputCS7.cleanAnotaton,'VE_CS7')==1 | strcmp(OutputCS7.cleanAnotaton,'Stalk_CS7')==1 );
indT4 = find(strcmp(OutputCS62.cleanAnotaton,'EmDisc_CS6')==1 | strcmp(OutputCS62.cleanAnotaton,'VE_CS6')==1 | strcmp(OutputCS62.cleanAnotaton,'EmDisc_Gast_CS6')==1 | strcmp(OutputCS62.cleanAnotaton,'EmDisc_Stalk_CS6')==1 | strcmp(OutputCS62.cleanAnotaton,'EmDisc_gast_CS6')==1 | strcmp(OutputCS62.cleanAnotaton,'Stalk_CS6')==1 | strcmp(OutputCS62.cleanAnotaton,'VE_CS6')==1);
m1 = min([OutputCS5.Ytrain(indT1),OutputCS6.Ytrain(indT2),OutputCS7.Ytrain(indT3)]);
m2 = max([OutputCS5.Ytrain(indT1),OutputCS6.Ytrain(indT2),OutputCS7.Ytrain(indT3)]);


subplot(2,2,1);
scatter3(OutputCS5c.cleanX(indT1,1),OutputCS5c.cleanX(indT1,2),OutputCS5c.cleanX(indT1,3),125,OutputCS5.Ytrain(indT1),'filled');
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(c1,d1)

subplot(2,2,2);
scatter3(OutputCS6c.cleanX(indT2,1),OutputCS6c.cleanX(indT2,2),OutputCS6c.cleanX(indT2,3),125,OutputCS6.Ytrain(indT2) , 'filled');
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(c2,d2)

subplot(2,2,3);
scatter3(OutputCS7c.cleanX(indT3,1),OutputCS7c.cleanX(indT3,2),OutputCS7c.cleanX(indT3,3),125,OutputCS7.Ytrain(indT3), 'filled' );
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(c3,d3)

subplot(2,2,4);
scatter3(OutputCS62c.cleanX(indT4,1),OutputCS62c.cleanX(indT4,2),OutputCS62c.cleanX(indT4,3),125,OutputCS62.Ytrain(indT4) ,'filled');
ax = gca(hh)
ax.CLim = [m1 m2];
axis equal
axis off
colormap(parula);
view(c2_2,d2_2)

print(hh,['Plots/Scatter_' gene{i} '_CS5-7.pdf'],'-dpdf','-r1000');

close all

    end
end