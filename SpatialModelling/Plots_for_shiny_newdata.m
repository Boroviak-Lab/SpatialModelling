addpath(genpath('./Functions'))

load('MissingSections/CS5_revisions_DIAGONAL_2.mat')
O1b = OBJ1_6;
load('./Data/Cross sections/CS6/CS6_revisions_CROSS_246.mat')
[O2b, sectcs6] = LoadCS6('246');
load('MissingSections/200319_CS7_build_full_stalkcut_FORDYLAN.mat')
[O3b, sectcs6] = LoadCS7('C2');

[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,Shots);
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS6');
[OutputCS6] = loadCS6Scaffold(D,Locations,Shots);
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS7');
[OutputCS7] = loadCS7Scaffold(D,Locations,Shots);

O1 = transformOBJ(O1b,pi,[0,0,1]);
O1 = transformOBJ(O1,'flipX');
O2 = transformOBJ(O2b,0.81*pi,[0,0,1]);
O2 = transformOBJ(O2,'flipX');
O3 = transformOBJ(O3b,-0.55*pi,[0,0,1]);


%%____________________________________________
%Offsets for moving around the Preimplantation
DeltaX = -100;
DeltaY = 300;

%CS5
O1.vertices(:,3) = O1.vertices(:,3) -180 ;

%CS6
O2.vertices(:,1) = O2.vertices(:,1) + 1500;
O2.vertices(:,3) = O2.vertices(:,3) - 180;

%CS7
O3.vertices(:,1) = -O3.vertices(:,1);
O3.vertices(:,1) = O3.vertices(:,1) + 4000;


load('MissingSections/lBlHR.mat')
O4.vertices = Quaternion3(pi*1.28,[0,1,0],O4.vertices); 
O4.vertices = O4.vertices*33.2;
O4.vertices(:,1) = O4.vertices(:,1) + 1800 + DeltaX;
O4.vertices(:,2) = O4.vertices(:,2) + 1000 + DeltaY;
O4.vertices(:,3) = O4.vertices(:,3) - 200;
f1 = patch('Faces',O4.objects(4).data.vertices,'Vertices',O4.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',O4.objects(12).data.vertices,'Vertices',O4.vertices,'FaceColor',[230,181,0]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O4.objects(8).data.vertices,'Vertices',O4.vertices,'FaceColor',[0,191,191]/255,'LineStyle','none','FaceAlpha',1);

load('MissingSections/eBlHR2.mat')
O5.vertices = Quaternion3(pi*0.95,[0,1,0],O5.vertices); 
O5.vertices = O5.vertices*33.2;
O5.vertices(:,1) = O5.vertices(:,1) + 1300 + DeltaX;
O5.vertices(:,2) = O5.vertices(:,2) + 950 + DeltaY;
O5.vertices(:,3) = O5.vertices(:,3) - 200;
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceColor',[0,230,230]/255,'LineStyle','none','FaceAlpha',1);

load('MissingSections/cMHR.mat')
O6.vertices = O6.vertices*33.2;
O6.vertices(:,1) = O6.vertices(:,1) + 800 + DeltaX;
O6.vertices(:,2) = O6.vertices(:,2) + 1020 + DeltaY;
%f1 = patch('Faces',O6.objects(4).data.vertices,'Vertices',O6.vertices,'FaceColor',[0,153,38]/255,'LineStyle','none','FaceAlpha',1);

load('MissingSections/ECHR.mat')
O7.vertices = O7.vertices*33.2;
O7.vertices(:,1) = O7.vertices(:,1) + 430 + DeltaX;
O7.vertices(:,2) = O7.vertices(:,2) + 1020 + DeltaY;

load('MissingSections/FCHR.mat')
O8.vertices = O8.vertices*33.2;
O8.vertices(:,1) = O8.vertices(:,1) + 90 + DeltaX;
O8.vertices(:,2) = O8.vertices(:,2) + 1000 + DeltaY;
O8.vertices(:,3) = O8.vertices(:,3) - 40;

load('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/LowRes/ZHR.mat')
O9.vertices = O9.vertices*33.2;
O9.vertices(:,1) = O9.vertices(:,1) - 300 + DeltaX;
O9.vertices(:,2) = O9.vertices(:,2) + 1000 + DeltaY;
O9.vertices(:,3) = O9.vertices(:,3) - 40;

%Now generate black background version (same as above)
close all
h = figure('visible', 'off');
DeltaX = -100;
DeltaY = 300;
%Amnion / 1 =Am/PGC
f1 = patch('Faces',O1.objects(4).data.vertices,'Vertices',  O1.vertices,'FaceColor', [135,123,214]/255,'LineStyle','none','FaceAlpha',1);
%EmDisc / 2 = EmDisc
f1 = patch('Faces',O1.objects(20).data.vertices,'Vertices',O1.vertices, 'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',1);
%VE / 4 = VE
f2 = patch('Faces',O1.objects(8).data.vertices,'Vertices',  O1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',1);
%SYS / 5 = SYS
f3 = patch('Faces',O1.objects(12).data.vertices,'Vertices',  O1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',1);
%ExMes / 6 = ExMes
f4 = patch('Faces',O1.objects(16).data.vertices,'Vertices',  O1.vertices,'FaceColor',[230,200,0]/255,'LineStyle','none','FaceAlpha',1);
%Tb / 7 = Tb
f5 = patch('Faces',O1.objects(24).data.vertices,'Vertices',  O1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',1);
%CS6
%Am / 
f1 = patch('Faces',O2.objects(20).data.vertices,'Vertices',  O2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',1);
%EmDisc
f2 = patch('Faces',O2.objects(4).data.vertices,'Vertices',  O2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',1);
%Stalk
f3 = patch('Faces',O2.objects(24).data.vertices,'Vertices',  O2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',1);
%PGC
f4 = patch('Faces',O2.objects(8).data.vertices,'Vertices',  O2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',1);
%VE
f5 = patch('Faces',O2.objects(16).data.vertices,'Vertices',  O2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',1);
%SYS
f6 = patch('Faces',O2.objects(12).data.vertices,'Vertices',  O2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',1);
%Tb
f7 = patch('Faces',O2.objects(28).data.vertices,'Vertices',  O2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',1);
%ExMes
f8 = patch('Faces',O2.objects(32).data.vertices,'Vertices',  O2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',1);


%Am
f1 = patch('Faces',O3.objects(4).data.vertices,'Vertices',  O3.vertices,'FaceColor',[113,8,166]/255,'LineStyle','none','FaceAlpha',1);
%EmDisc
f2 = patch('Faces',O3.objects(8).data.vertices,'Vertices',  O3.vertices,'FaceColor',[191,113,4]/255,'LineStyle','none','FaceAlpha',1);
%Stalk
f3 = patch('Faces',O3.objects(20).data.vertices,'Vertices',  O3.vertices,'FaceColor',[2,51,191]/255,'LineStyle','none','FaceAlpha',1);
%VE
f4 = patch('Faces',O3.objects(12).data.vertices,'Vertices',  O3.vertices,'FaceColor',[191,60,4]/255,'LineStyle','none','FaceAlpha',1);
%SYS
f5 = patch('Faces',O3.objects(16).data.vertices,'Vertices',  O3.vertices,'FaceColor',[26,8,115]/255,'LineStyle','none','FaceAlpha',1);
%Tb
f6 = patch('Faces',O3.objects(30).data.vertices,'Vertices',  O3.vertices,'FaceColor',[96, 56, 19]/255,'LineStyle','none','FaceAlpha',1);
%ExMes
f7 = patch('Faces',O3.objects(34).data.vertices,'Vertices',  O3.vertices,'FaceColor',[150,119,0]/255,'LineStyle','none','FaceAlpha',1);
%Pre
f1 = patch('Faces',O4.objects(4).data.vertices,'Vertices',O4.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',O4.objects(12).data.vertices,'Vertices',O4.vertices,'FaceColor',[230,181,0]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O4.objects(8).data.vertices,'Vertices',O4.vertices,'FaceColor',[0,191,191]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceColor',[0,230,230]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O6.objects(4).data.vertices,'Vertices',O6.vertices,'FaceColor',[0,153,38]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O7.objects(10).data.vertices,'Vertices',O7.vertices,'FaceColor',[0,191,48]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O8.objects(9).data.vertices,'Vertices',O8.vertices,'FaceColor',[72,191,0]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O9.objects(4).data.vertices,'Vertices',O9.vertices,'FaceColor',[86,230,0]/255,'LineStyle','none','FaceAlpha',1);
hold on
%Scale bars
plot(linspace(4600,4600+400,100),-400*ones(1,100),'w','LineWidth',5)
xlim([-850,6000])
ylim([-800,1900])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./All/Lineage_bBG_OFF.png'],'-r800')

text(4710,-340,{'200 um'},'FontSize',8,'Color',[1, 1 ,1])
 text(-600+ DeltaX,1300+ DeltaY,{'Zygote'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-310,-360,100)+ DeltaX,linspace(1016,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(-120+ DeltaX,1300+ DeltaY,{'4-cell'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(65,15,100)+ DeltaX,linspace(1030,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(310+ DeltaX,1300+ DeltaY,{'8-cell'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(472,500,100)+ DeltaX,linspace(1065,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(750+ DeltaX,1300+ DeltaY,{'Morula'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(846,906,100)+ DeltaX,linspace(1095,1220,100)+ DeltaY,'w-','LineWidth',1)
 text(1200+ DeltaX,1300+ DeltaY,{'ICM'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1326,1325,100)+ DeltaX,linspace(1077,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(1500+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1320,1480,100)+ DeltaX,linspace(915,850,100)+ DeltaY,'w-','LineWidth',1)
 text(2000+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1905,2000,100)+ DeltaX,linspace(882,850,100)+ DeltaY,'w-','LineWidth',1)
 text(1700+ DeltaX,1300+ DeltaY,{'Epiblast'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1806,1790,100)+ DeltaX,linspace(1124,1240,100)+ DeltaY,'w-','LineWidth',1)
 text(2100+ DeltaX,1100+ DeltaY,{'Hypoblast'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1932,2090,100)+ DeltaX,linspace(1067,1065,100)+ DeltaY,'w-','LineWidth',1)
 
 text(-680,150,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-206,-410,100),linspace(85,130,100),'w-','LineWidth',1)
 text(-230,275,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-74,-180,100),linspace(110,220,100),'w-','LineWidth',1)
 text(0,250,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(17,0,100),linspace(190,70,100),'w-','LineWidth',1)
 text(280,-240,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(44,250,100),linspace(-23,-200,100),'w-','LineWidth',1)
 text(-480,-330,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-233,-380,100),linspace(-40,-220,100),'w-','LineWidth',1)
 text(-110,-360,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(0,-10,100),linspace(-33,-310,100),'w-','LineWidth',1)
%CS6
text(900,-340,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1135,1000,100),linspace(-220,-300,100),'w-','LineWidth',1)
text(2000,-300,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1839,2100,100),linspace(-76,-250,100),'w-','LineWidth',1)
text(1940,50,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1829,1930,100),linspace(-10,50,100),'w-','LineWidth',1)
text(1700,-400,{'PGC'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1770,1820,100),linspace(-70,-350,100),'w-','LineWidth',1)
text(800,200,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1200,1000,100),linspace(120,180,100),'w-','LineWidth',1)
text(1100,400,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1278,1300,100),linspace(150,350,100),'w-','LineWidth',1)
text(1600,400,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1572,1600,100),linspace(240,350,100),'w-','LineWidth',1)
text(1820,220,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1657,1820,100),linspace(91,220,100),'w-','LineWidth',1)

text(2900,-350,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(3147,3000,100),linspace(-143,-300,100),'w-','LineWidth',1)
text(5300,600,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(5152,5300,100),linspace(557,600,100),'w-','LineWidth',1)
text(4200,1500,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4234,4300,100),linspace(1147,1450,100),'w-','LineWidth',1)
text(2370,480,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(2690,2500,100),linspace(577,510,100),'w-','LineWidth',1)
plot(linspace(4843,5400,100),linspace(104,200,100),'w-','LineWidth',1)
text(5400,200,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
text(5360,400,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4673,5304,100),linspace(160,390,100),'w-','LineWidth',1)


 text(-1115,900,{'Carnegie Stage 1'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
 text(83,900,{'Carnegie Stage 2'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
 text(1280,900,{'Carnegie Stage 3'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
 text(-400,-550,{'Carnegie Stage 5'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
 text(1100,-550,{'Carnegie Stage 6'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
 text(3200,-550,{'Carnegie Stage 7'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
 
xlim([-850,6000])
ylim([-800,1900])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./All/Lineage_bBG_ON.png'],'-r800')

close all

%Now the CS5 only 
h = figure('visible', 'off');
f1 = patch('Faces',O1.objects(4).data.vertices,'Vertices',  O1.vertices,'FaceColor', [135,123,214]/255,'LineStyle','none','FaceAlpha',1);
%EmDisc / 2 = EmDisc
f1 = patch('Faces',O1.objects(20).data.vertices,'Vertices',O1.vertices, 'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',1);
%VE / 4 = VE
f2 = patch('Faces',O1.objects(8).data.vertices,'Vertices',  O1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',1);
%SYS / 5 = SYS
f3 = patch('Faces',O1.objects(12).data.vertices,'Vertices',  O1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',1);
%ExMes / 6 = ExMes
f4 = patch('Faces',O1.objects(16).data.vertices,'Vertices',  O1.vertices,'FaceColor',[230,200,0]/255,'LineStyle','none','FaceAlpha',1);
%Tb / 7 = Tb
f5 = patch('Faces',O1.objects(24).data.vertices,'Vertices',  O1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O1.objects(4).data.vertices,'Vertices',  O1.vertices,'FaceColor', [135,123,214]/255,'LineStyle','none','FaceAlpha',1);
%EmDisc / 2 = EmDisc
f1 = patch('Faces',O1.objects(20).data.vertices,'Vertices',O1.vertices, 'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',1);
%VE / 4 = VE
f2 = patch('Faces',O1.objects(8).data.vertices,'Vertices',  O1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',1);
%SYS / 5 = SYS
f3 = patch('Faces',O1.objects(12).data.vertices,'Vertices',  O1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',1);
%ExMes / 6 = ExMes
f4 = patch('Faces',O1.objects(16).data.vertices,'Vertices',  O1.vertices,'FaceColor',[230,200,0]/255,'LineStyle','none','FaceAlpha',1);
%Tb / 7 = Tb
f5 = patch('Faces',O1.objects(24).data.vertices,'Vertices',  O1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',1);
hold on
xlim([-700 600])
ylim([-750 550])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS5/Lineage_CS5_bBG_OFF.png'],'-r800')
text(-530,150,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-206,-410,100),linspace(110,130,100),'w-','LineWidth',1)
 text(-230,265,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-74,-180,100),linspace(110,220,100),'w-','LineWidth',1)
 text(0,240,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(17,0,100),linspace(190,70,100),'w-','LineWidth',1)
 text(280,-240,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(44,250,100),linspace(-23,-200,100),'w-','LineWidth',1)
 text(-460,-290,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-233,-380,100),linspace(-40,-220,100),'w-','LineWidth',1)
 text(-100,-360,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(0,-10,100),linspace(-33,-310,100),'w-','LineWidth',1)
text(-200,490,{'LINEAGE'},'FontSize',22,'Color',[1, 1 ,1])
text(-300,410,{'Carnegie Stage 5'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS5/Lineage_CS5_bBG_ON.png'],'-r800')

close all
%CS6
h = figure('visible', 'off');
f1 = patch('Faces',O2.objects(20).data.vertices,'Vertices',  O2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',1);
%EmDisc
f2 = patch('Faces',O2.objects(4).data.vertices,'Vertices',  O2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',1);
%Stalk
f3 = patch('Faces',O2.objects(24).data.vertices,'Vertices',  O2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',1);
%PGC
f4 = patch('Faces',O2.objects(8).data.vertices,'Vertices',  O2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',1);
%VE
f5 = patch('Faces',O2.objects(16).data.vertices,'Vertices',  O2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',1);
%SYS
f6 = patch('Faces',O2.objects(12).data.vertices,'Vertices',  O2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',1);
%Tb
f7 = patch('Faces',O2.objects(28).data.vertices,'Vertices',  O2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',1);
%ExMes
f8 = patch('Faces',O2.objects(32).data.vertices,'Vertices',  O2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',1);
hold on
xlim([780 2400])
ylim([-870 750])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS6/Lineage_CS6_bBG_OFF.png'],'-r800')
text(900,-340,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1135,1000,100),linspace(-220,-300,100),'w-','LineWidth',1)
text(2040,-300,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1839,2100,100),linspace(-76,-250,100),'w-','LineWidth',1)
text(1940,50,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1829,1930,100),linspace(-10,50,100),'w-','LineWidth',1)
text(1740,-400,{'PGC'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1770,1820,100),linspace(-70,-350,100),'w-','LineWidth',1)
text(800,200,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1200,1000,100),linspace(120,180,100),'w-','LineWidth',1)
text(1150,400,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1278,1300,100),linspace(150,350,100),'w-','LineWidth',1)
text(1570,400,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1572,1600,100),linspace(240,350,100),'w-','LineWidth',1)
text(1820,220,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1657,1820,100),linspace(91,220,100),'w-','LineWidth',1)

text(1300,675,{'LINEAGE'},'FontSize',22,'Color',[1, 1 ,1])
text(1200,575,{'Carnegie Stage 6'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS6/Lineage_CS6_bBG_ON.png'],'-r800')

close all


%CS7
h = figure('visible', 'off');
%Am
f1 = patch('Faces',O3.objects(4).data.vertices,'Vertices',  O3.vertices,'FaceColor',[113,8,166]/255,'LineStyle','none','FaceAlpha',1);
%EmDisc
f2 = patch('Faces',O3.objects(8).data.vertices,'Vertices',  O3.vertices,'FaceColor',[191,113,4]/255,'LineStyle','none','FaceAlpha',1);
%Stalk
f3 = patch('Faces',O3.objects(20).data.vertices,'Vertices',  O3.vertices,'FaceColor',[2,51,191]/255,'LineStyle','none','FaceAlpha',1);
%VE
f4 = patch('Faces',O3.objects(12).data.vertices,'Vertices',  O3.vertices,'FaceColor',[191,60,4]/255,'LineStyle','none','FaceAlpha',1);
%SYS
f5 = patch('Faces',O3.objects(16).data.vertices,'Vertices',  O3.vertices,'FaceColor',[26,8,115]/255,'LineStyle','none','FaceAlpha',1);
%Tb
f6 = patch('Faces',O3.objects(30).data.vertices,'Vertices',  O3.vertices,'FaceColor',[96, 56, 19]/255,'LineStyle','none','FaceAlpha',1);
%ExMes
f7 = patch('Faces',O3.objects(34).data.vertices,'Vertices',  O3.vertices,'FaceColor',[150,119,0]/255,'LineStyle','none','FaceAlpha',1);
hold on
xlim([2220 6000])
ylim([-1180 2600])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS7/Lineage_CS7_bBG_OFF.png'],'-r1000')
%print('-dpng',['~/Desktop/Lineage_CS7_bBG_OFF.png'],'-r800')
text(2900,-350,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
hold on
plot(linspace(3147,3000,100),linspace(-143,-300,100),'w-','LineWidth',1)
text(5300,600,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(5152,5300,100),linspace(557,600,100),'w-','LineWidth',1)
text(4200,1500,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4234,4300,100),linspace(1147,1450,100),'w-','LineWidth',1)
text(2370,480,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(2690,2500,100),linspace(577,510,100),'w-','LineWidth',1)
plot(linspace(4843,5400,100),linspace(104,200,100),'w-','LineWidth',1)
text(5400,200,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
text(5360,400,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4673,5304,100),linspace(160,390,100),'w-','LineWidth',1)

text(3700,2430,{'LINEAGE'},'FontSize',22,'Color',[1, 1 ,1])
text(3400,2220,{'Carnegie Stage 7'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS7/Lineage_CS7_bBG_ON.png'],'-r1000')

close all

%Pre: CS3
h = figure('visible', 'off');

f1 = patch('Faces',O4.objects(4).data.vertices,'Vertices',O4.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',O4.objects(12).data.vertices,'Vertices',O4.vertices,'FaceColor',[230,181,0]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O4.objects(8).data.vertices,'Vertices',O4.vertices,'FaceColor',[0,191,191]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceColor',[0,230,230]/255,'LineStyle','none','FaceAlpha',1);
xlim([900,2300])
ylim([700,2100])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS3/Lineage_CS3_OFF.png'],'-r1000')
hold on
text(1300+ DeltaX,1300+ DeltaY,{'ICM'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1326,1325,100)+ DeltaX,linspace(1077,1230,100)+ DeltaY,'w-','LineWidth',1)
text(1500+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1320,1480,100)+ DeltaX,linspace(915,850,100)+ DeltaY,'w-','LineWidth',1)
text(2000+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1905,2000,100)+ DeltaX,linspace(882,850,100)+ DeltaY,'w-','LineWidth',1)
text(1700+ DeltaX,1300+ DeltaY,{'Epiblast'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1806,1790,100)+ DeltaX,linspace(1124,1240,100)+ DeltaY,'w-','LineWidth',1)
text(2100+ DeltaX,1050+ DeltaY,{'Hypoblast'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1932,2090,100)+ DeltaX,linspace(1067,1065,100)+ DeltaY,'w-','LineWidth',1)
text(1400,2050,{'LINEAGE'},'FontSize',22,'Color',[1, 1 ,1])
text(1300,1960,{'Carnegie Stage 3'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS3/Lineage_CS3_ON.png'],'-r800')

close all

%Pre CS1-2
h = figure('visible', 'off');

f1 = patch('Faces',O6.objects(4).data.vertices,'Vertices',O6.vertices,'FaceColor',[0,153,38]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O7.objects(10).data.vertices,'Vertices',O7.vertices,'FaceColor',[0,191,48]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O8.objects(9).data.vertices,'Vertices',O8.vertices,'FaceColor',[72,191,0]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O9.objects(4).data.vertices,'Vertices',O9.vertices,'FaceColor',[86,230,0]/255,'LineStyle','none','FaceAlpha',1);
xlim([-600,1200])
ylim([200,2000])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS1-2/Lineage_CS1-2_bBG_OFF.png'],'-r800')
hold on
text(-460+ DeltaX,1300+ DeltaY,{'Zygote'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(-310,-360,100)+ DeltaX,linspace(1016,1230,100)+ DeltaY,'w-','LineWidth',1)
text(-20+ DeltaX,1300+ DeltaY,{'4-cell'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(65,15,100)+ DeltaX,linspace(1030,1230,100)+ DeltaY,'w-','LineWidth',1)
text(350+ DeltaX,1300+ DeltaY,{'8-cell'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(472,500,100)+ DeltaX,linspace(1065,1230,100)+ DeltaY,'w-','LineWidth',1)
text(700+ DeltaX,1300+ DeltaY,{'Morula'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(846,906,100)+ DeltaX,linspace(1095,1220,100)+ DeltaY,'w-','LineWidth',1)
text(120,1950,{'LINEAGE'},'FontSize',22,'Color',[1, 1 ,1])
text(-50,1850,{'Carnegie Stages 1-2'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['~/Desktop/OnlineResource/CS1-2/Lineage_CS1-2_bBG_ON.png'],'-r800')


close all

%___________________________________________
%Now do some GP interpolations for genes

IDs= importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/SpatialModelling/Data/IDs_and_locations.csv')
Dexp = importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/SpatialModelling/Data/X.csv')


%list = {'SOX17','SOX2','POU5F1','T','TFAP2C','NANOG','PRDM14','LEF1','ID2','BMP4','TDGF1','NODAL','WNT3','WNT8A','TFAP2A','MIXL1','CER1','DPPA5','FBXO2','PDGFRA','PDGFA','KDR','VEGFA','DAZL','KLF4','PRAME','MAEL','WNT5A','WNT5B','WNT6','TCF4','ID1','ID3','CHRD','FST','MIXL1','GATA6','CDH2','LEFTY2','LHX1','OTX2','HHEX','SDCC4','BMP2','APOA1','GC','IHH','ARL13B','PTCH2','SMO','GPC4','GLI1'}

Lind1 = find(strcmp(IDs.textdata(2:end,4),'Zy_CS1')==1 );
Lind2 = find(strcmp(IDs.textdata(2:end,4),'4-cell_CS2')==1  );
Lind3 = find(strcmp(IDs.textdata(2:end,4),'8-cell_CS2')==1  );
Lind4 = find(strcmp(IDs.textdata(2:end,4),'cMor_CS3')==1  );
Lind5 = find(strcmp(IDs.textdata(2:end,4),'ICM_CS3')==1  );
Lind6 = find(strcmp(IDs.textdata(2:end,4),'Epi_CS3')==1  );
Lind7 = find(strcmp(IDs.textdata(2:end,4),'Tb_CS3')==1 );
Lind8 = find(strcmp(IDs.textdata(2:end,4),'Hyp_CS3')==1  );

load(['~/Desktop/Hyp1_batch1.mat']);
hcs6=load(['~/Desktop/CS6_Hyp1_batch1.mat']);
Hyp3 = hcs6.Hyp1;
hcs7=load(['~/Desktop/CS7_Hyp1_batch1.mat']);
Hyp5 = hcs7.Hyp1;
list = find(strcmp(Dexp.textdata(2:end,1),'SOX2')==1);
genelist = Dexp.textdata(2:end,1);

for i = 1:length(genelist)    
    try
        II(i,1) = strfind(genelist{i},'ENSCJAG'); 
    catch
        II(i,1) = 0;         
    end    
end

genelist2 = genelist(find(II==0));
%genelist2 = {'POU5F1','NANOG','SOX2','SFRP2','DNMT3B','T','MIXL1','EOMES','TFAP2C','TFAP2A','VTCN1','SOX17','PRDM1','NANOS3','PRDM14','DND1','FOXA2'};

[OutputCS5] = MarmosetGP_CS5Opt(D,OutputCS5,1);
[OutputCS6] = MarmosetGP_CS6_v3Opt(D,OutputCS6,1);
[OutputCS7] = MarmosetGP_CS7_v3Opt(D,OutputCS7,1);


for i = 8000:length(genelist2)
close all    


if(int64(i/1000)==(i/1000))
    disp(['Step ' num2str(i)])
end
    try            
list = find(strcmp(Dexp.textdata(2:end,1),genelist2{i})==1);

            
%Generate GP models using preprocessed hyperparams            
[OutputCS5] = MarmosetGPInfer_CS5Opt(OutputCS5,O1b,list,D,Hyp1(list,:),Hyp1(list,:));
[OutputCS6] = MarmosetGPInfer_CS6_v3Opt(OutputCS6,O2b,list,D,Hyp3(list,:),Hyp3(list,:));
[OutputCS7] = MarmosetGPInfer_CS7_v3Opt(OutputCS7,O3b,list,D,Hyp5(list,:),Hyp5(list,:));       
            
%[Output1] = Marmoset3D_CS5_surface(DexpCS5,OutputCS5,'Base2',list{i});
%[Output2] = Marmoset3D_CS6_surface(DexpCS6,OutputCS6,'Base2',list{i});
%[Output3] = Marmoset3D_CS7_surface(DexpCS7,OutputCS7,'Base2',list{i});

minV = min([min(OutputCS5.cLim),min(OutputCS6.cLim),min(OutputCS7.cLim)]);
maxV = max([max(OutputCS5.cLim),max(OutputCS6.cLim),max(OutputCS7.cLim)]);
%maxV = max([Output1.m_0;Output1.m_1;Output1.m_2;Output1.m_3;Output1.m_4;Output1.m_5;   Output2.m_0;Output2.m_1;Output2.m_2; Output2.m_3; Output2.m_4;Output2.m_5;Output2.m_6;    Output3.m_0;Output3.m_1;Output3.m_2;Output3.m_3;Output3.m_5]);
%minV = min([Output1.m_0;Output1.m_1;Output1.m_2;Output1.m_3;Output1.m_4;Output1.m_5;   Output2.m_0;Output2.m_1;Output2.m_2; Output2.m_3; Output2.m_4;Output2.m_5;Output2.m_6;    Output3.m_0;Output3.m_1;Output3.m_2;Output3.m_3;Output3.m_5]);

M1 = mean(Dexp.data(list,Lind1));
M2 = mean(Dexp.data(list,Lind2));
M3 = mean(Dexp.data(list,Lind3));
M4 = mean(Dexp.data(list,Lind4));
M5 = mean(Dexp.data(list,Lind5));
M6 = mean(Dexp.data(list,Lind6));
M7 = mean(Dexp.data(list,Lind7));
M8 = mean(Dexp.data(list,Lind8));

close all
h = figure('visible', 'off');
%Amnion / 1 =Am/PGC
f1 = patch('Faces',O1.objects(4).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
%EmDisc / 2 = EmDisc
f1 = patch('Faces',O1.objects(20).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
%VE / 4 = VE
f2 = patch('Faces',O1.objects(8).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
%SYS / 5 = SYS
f3 = patch('Faces',O1.objects(12).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
%ExMes / 6 = ExMes
f4 = patch('Faces',O1.objects(16).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
%Tb / 7 = Tb
f5 = patch('Faces',O1.objects(24).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
%CS6
%Am / 
f1 = patch('Faces',O2.objects(20).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
%EmDisc
f2 = patch('Faces',O2.objects(4).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
%Stalk
f3 = patch('Faces',O2.objects(24).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m8,'FaceColor','interp','LineStyle','none');
%PGC
f4 = patch('Faces',O2.objects(8).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
%VE
f5 = patch('Faces',O2.objects(16).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
%SYS
f6 = patch('Faces',O2.objects(12).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
%Tb
f7 = patch('Faces',O2.objects(28).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
%ExMes
f8 = patch('Faces',O2.objects(32).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');

%CS7
%Tb
f1 = patch('Faces',O3.objects(4).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m7,'FaceColor','interp','LineStyle','none');
%SYS
f2 = patch('Faces',O3.objects(8).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m5,'FaceColor','interp','LineStyle','none');
%EmDisc
f3 = patch('Faces',O3.objects(20).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m2,'FaceColor','interp','LineStyle','none');
%VE
%f4 = patch('Faces',O3.objects(12).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m1,'FaceColor','interp','LineStyle','none');
%Am
f5 = patch('Faces',O3.objects(16).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m1,'FaceColor','interp','LineStyle','none');
%Stalk
f6 = patch('Faces',O3.objects(30).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m8,'FaceColor','interp','LineStyle','none');
%ExMes
f7 = patch('Faces',O3.objects(34).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m6,'FaceColor','interp','LineStyle','none');

%Preimplantation

f1 = patch('Faces',O4.objects(4).data.vertices,'Vertices',O4.vertices,'FaceVertexCData',M7*ones(length(O4.vertices),1),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',O4.objects(12).data.vertices,'Vertices',O4.vertices,'FaceVertexCData',M8*ones(length(O4.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O4.objects(8).data.vertices,'Vertices',O4.vertices,'FaceVertexCData',M6*ones(length(O4.vertices),1),'FaceColor','interp','LineStyle','none');

f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',M7*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',M7*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',M5*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
%f1 = patch('Faces',O4.objects(4).data.vertices,'Vertices',O4.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
%f3 = patch('Faces',O4.objects(12).data.vertices,'Vertices',O4.vertices,'FaceColor',[230,181,0]/255,'LineStyle','none','FaceAlpha',1);
%f2 = patch('Faces',O4.objects(8).data.vertices,'Vertices',O4.vertices,'FaceColor',[0,191,191]/255,'LineStyle','none','FaceAlpha',1);
%f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceColor',[191,4,137]/255,'LineStyle','none','FaceAlpha',1);
%f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceColor',[0,230,230]/255,'LineStyle','none','FaceAlpha',1);
f1 = patch('Faces',O6.objects(4).data.vertices,'Vertices',O6.vertices,'FaceVertexCData',M4*ones(length(O6.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O7.objects(10).data.vertices,'Vertices',O7.vertices,'FaceVertexCData',M3*ones(length(O7.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O8.objects(9).data.vertices,'Vertices',O8.vertices,'FaceVertexCData',M2*ones(length(O8.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O9.objects(4).data.vertices,'Vertices',O9.vertices,'FaceVertexCData',M1*ones(length(O9.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
%Add scale bars
%plot(linspace(1500,2000,100)+ DeltaX,600*ones(1,100)+ DeltaY,'w','LineWidth',5)
%text(1655+DeltaX,670+DeltaY,{'100 um'},'FontSize',8,'Color',[1, 1 ,1])
%plot(linspace(4000,5000,100),-400*ones(1,100),'w','LineWidth',5)
%text(4710,-340,{'100 um'},'FontSize',8,'Color',[1, 1 ,1])
%Add labels
plot(linspace(4600,4600+400,100),-400*ones(1,100),'w','LineWidth',5) 
xlim([-850,6000])
ylim([-800,1900])
%Axis properties
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf,'Color','k')
set(gca,'clim',[minV maxV] )
cb=colorbar;
cb.Position = cb.Position + [0 0.5 -0.0186 -0.6]
cb.Color = 'w';
cb.FontSize = 12;
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(2),cy(end)])
end
set(gca,'fontsize', 12)
set(gcf, 'InvertHardcopy', 'off')
axis off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./All/' genelist{list} '_bBG_OFF.png'],'-r800')

text(4710,-340,{'200 um'},'FontSize',8,'Color',[1, 1 ,1])


h0 = text(5350,850,{'mRNA level'},'fontweight','FontSize',10,'Color',[1, 1 ,1])
set(h0,'Rotation',90);

text(4710,-340,{'200 um'},'FontSize',8,'Color',[1, 1 ,1])
 text(-600+ DeltaX,1300+ DeltaY,{'Zygote'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-310,-360,100)+ DeltaX,linspace(1016,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(-120+ DeltaX,1300+ DeltaY,{'4-cell'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(65,15,100)+ DeltaX,linspace(1030,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(310+ DeltaX,1300+ DeltaY,{'8-cell'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(472,500,100)+ DeltaX,linspace(1065,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(750+ DeltaX,1300+ DeltaY,{'Morula'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(846,906,100)+ DeltaX,linspace(1095,1220,100)+ DeltaY,'w-','LineWidth',1)
 text(1200+ DeltaX,1300+ DeltaY,{'ICM'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1326,1325,100)+ DeltaX,linspace(1077,1230,100)+ DeltaY,'w-','LineWidth',1)
 text(1500+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1320,1480,100)+ DeltaX,linspace(915,850,100)+ DeltaY,'w-','LineWidth',1)
 text(2000+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1905,2000,100)+ DeltaX,linspace(882,850,100)+ DeltaY,'w-','LineWidth',1)
 text(1700+ DeltaX,1300+ DeltaY,{'Epiblast'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1806,1790,100)+ DeltaX,linspace(1124,1240,100)+ DeltaY,'w-','LineWidth',1)
 text(2100+ DeltaX,1100+ DeltaY,{'Hypoblast'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(1932,2090,100)+ DeltaX,linspace(1067,1065,100)+ DeltaY,'w-','LineWidth',1)
 text(-680,150,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-206,-410,100),linspace(85,130,100),'w-','LineWidth',1)
 text(-230,275,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-74,-180,100),linspace(110,220,100),'w-','LineWidth',1)
 text(0,250,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(17,0,100),linspace(190,70,100),'w-','LineWidth',1)
 text(280,-240,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(44,250,100),linspace(-23,-200,100),'w-','LineWidth',1)
 text(-480,-330,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(-233,-380,100),linspace(-40,-220,100),'w-','LineWidth',1)
 text(-110,-360,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
 plot(linspace(0,-10,100),linspace(-33,-310,100),'w-','LineWidth',1)
%CS6
text(900,-340,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1135,1000,100),linspace(-220,-300,100),'w-','LineWidth',1)
text(2000,-300,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1839,2100,100),linspace(-76,-250,100),'w-','LineWidth',1)
text(1940,50,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1829,1930,100),linspace(-10,50,100),'w-','LineWidth',1)
text(1700,-400,{'PGC'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1770,1820,100),linspace(-70,-350,100),'w-','LineWidth',1)
text(800,200,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1200,1000,100),linspace(120,180,100),'w-','LineWidth',1)
text(1100,400,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1278,1300,100),linspace(150,350,100),'w-','LineWidth',1)
text(1600,400,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1572,1600,100),linspace(240,350,100),'w-','LineWidth',1)
text(1820,220,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1657,1820,100),linspace(91,220,100),'w-','LineWidth',1)
text(2900,-350,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(3147,3000,100),linspace(-143,-300,100),'w-','LineWidth',1)
text(5300,600,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(5152,5300,100),linspace(557,600,100),'w-','LineWidth',1)
text(4200,1500,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4234,4300,100),linspace(1147,1450,100),'w-','LineWidth',1)
text(2370,480,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(2690,2500,100),linspace(577,510,100),'w-','LineWidth',1)
plot(linspace(4843,5400,100),linspace(104,200,100),'w-','LineWidth',1)
text(5400,200,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
text(5360,400,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4673,5304,100),linspace(160,390,100),'w-','LineWidth',1)
text(-1115,900,{'Carnegie Stage 1'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
text(83,900,{'Carnegie Stage 2'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
text(1280,900,{'Carnegie Stage 3'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
text(-400,-550,{'Carnegie Stage 5'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
text(1100,-550,{'Carnegie Stage 6'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
text(3200,-550,{'Carnegie Stage 7'},'fontweight','bold','FontSize',12,'Color',[1, 1 ,1])
print('-dpng',['./All/' genelist{list} '_bBG_ON.png'],'-r800')
close all

%Now the CS5 only 
h = figure('visible', 'off');
f1 = patch('Faces',O1.objects(4).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m1,'FaceColor','interp','LineStyle','none');
%EmDisc / 2 = EmDisc
f1 = patch('Faces',O1.objects(20).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m2,'FaceColor','interp','LineStyle','none');
%VE / 4 = VE
f2 = patch('Faces',O1.objects(8).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m4,'FaceColor','interp','LineStyle','none');
%SYS / 5 = SYS
f3 = patch('Faces',O1.objects(12).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m5,'FaceColor','interp','LineStyle','none');
%ExMes / 6 = ExMes
f4 = patch('Faces',O1.objects(16).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m6,'FaceColor','interp','LineStyle','none');
%Tb / 7 = Tb
f5 = patch('Faces',O1.objects(24).data.vertices,'Vertices',  O1.vertices,'FaceVertexCData',OutputCS5.m7,'FaceColor','interp','LineStyle','none');
hold on
xlim([-700 600])
ylim([-750 550])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
%axis off
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gca,'clim',[minV maxV] )
cb=colorbar;
cb.Color = 'w';
cb.FontSize = 16;
cb.Position = cb.Position + [-0.1 0.5 -0 -0.6]
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(2),cy(end)])
end
set(gcf, 'InvertHardcopy', 'off')
axis off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
print('-dpng',['./CS5/' genelist{list} '_CS5_bBG_OFF.png'],'-r800')


h0 = text(300,70,{'mRNA level'},'fontweight','FontSize',18,'Color',[1, 1 ,1])
set(h0,'Rotation',90);
text(-530,150,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
hold on
plot(linspace(-206,-410,100),linspace(110,130,100),'w-','LineWidth',1)
text(-230,265,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(-74,-180,100),linspace(110,220,100),'w-','LineWidth',1)
text(0,240,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(17,0,100),linspace(190,70,100),'w-','LineWidth',1)
text(280,-240,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(44,250,100),linspace(-23,-200,100),'w-','LineWidth',1)
text(-460,-290,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(-233,-380,100),linspace(-40,-220,100),'w-','LineWidth',1)
text(-100,-360,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(0,-10,100),linspace(-33,-310,100),'w-','LineWidth',1)
text(-200,490,genelist{list},'FontSize',22,'Color',[1, 1 ,1])
text(-300,410,{'Carnegie Stage 5'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS5/' genelist{list} '_CS5_bBG_ON.png'],'-r800')
close all

%CS6
h = figure('visible', 'off');
f1 = patch('Faces',O2.objects(20).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m1,'FaceColor','interp','LineStyle','none');
%EmDisc
f2 = patch('Faces',O2.objects(4).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m2,'FaceColor','interp','LineStyle','none');
%Stalk
f3 = patch('Faces',O2.objects(24).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m8,'FaceColor','interp','LineStyle','none');
%PGC
f4 = patch('Faces',O2.objects(8).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m3,'FaceColor','interp','LineStyle','none');
%VE
f5 = patch('Faces',O2.objects(16).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m4,'FaceColor','interp','LineStyle','none');
%SYS
f6 = patch('Faces',O2.objects(12).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m5,'FaceColor','interp','LineStyle','none');
%Tb
f7 = patch('Faces',O2.objects(28).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m7,'FaceColor','interp','LineStyle','none');
%ExMes
f8 = patch('Faces',O2.objects(32).data.vertices,'Vertices',  O2.vertices,'FaceVertexCData',OutputCS6.m6,'FaceColor','interp','LineStyle','none');
hold on
xlim([780 2400])
ylim([-870 750])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gca,'clim',[minV maxV] )
cb=colorbar;
cb.Color = 'w';
cb.FontSize = 16;
cb.Position = cb.Position + [-0.05 0.45 -0 -0.6]
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
set(cb,'YTick',[cy(1),cy(2),cy(end)])    
end
set(gcf, 'InvertHardcopy', 'off')
%Save
axis off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
print('-dpng',['./CS6/' genelist{list} '_CS6_bBG_OFF.png'],'-r800')
%print('-dpng',['~/Desktop/' list{i} '_CS6_bBG_OFF.png'],'-r800')
h0 = text(2100,70,{'mRNA level'},'fontweight','FontSize',18,'Color',[1, 1 ,1])
set(h0,'Rotation',90);
text(900,-340,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1135,1000,100),linspace(-220,-300,100),'w-','LineWidth',1)
text(2040,-300,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1839,2100,100),linspace(-76,-250,100),'w-','LineWidth',1)
text(1940,50,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1829,1930,100),linspace(-10,50,100),'w-','LineWidth',1)
text(1740,-400,{'PGC'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1770,1820,100),linspace(-70,-350,100),'w-','LineWidth',1)
text(800,200,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1200,1000,100),linspace(120,180,100),'w-','LineWidth',1)
text(1150,400,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1278,1300,100),linspace(150,350,100),'w-','LineWidth',1)
text(1570,400,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1572,1600,100),linspace(240,350,100),'w-','LineWidth',1)
text(1820,220,{'VE'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(1657,1820,100),linspace(91,220,100),'w-','LineWidth',1)
text(1300,675,genelist{list},'FontSize',22,'Color',[1, 1 ,1])
text(1200,575,{'Carnegie Stage 6'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS6/' genelist{list} '_CS6_bBG_ON.png'],'-r800')
close all

%Now run CS7 only
h = figure('visible', 'off');
f1 = patch('Faces',O3.objects(4).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m8,'FaceColor','interp','LineStyle','none');
%SYS
f2 = patch('Faces',O3.objects(8).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m5,'FaceColor','interp','LineStyle','none');
%EmDisc
f3 = patch('Faces',O3.objects(20).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m2,'FaceColor','interp','LineStyle','none');
%VE
%f4 = patch('Faces',O3.objects(12).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m1,'FaceColor','interp','LineStyle','none');
%Am
f5 = patch('Faces',O3.objects(16).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m1,'FaceColor','interp','LineStyle','none');
%Stalk
f6 = patch('Faces',O3.objects(30).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m8,'FaceColor','interp','LineStyle','none');
%ExMes
f7 = patch('Faces',O3.objects(34).data.vertices,'Vertices',  O3.vertices,'FaceVertexCData',OutputCS7.m6,'FaceColor','interp','LineStyle','none');
hold on
xlim([2220 6000])
ylim([-1180 2600])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
%axis off
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gca,'clim',[minV maxV] )
cb=colorbar;
cb.Color = 'w';
cb.FontSize = 16;
cb.Position = cb.Position + [-0.05 0.45 -0 -0.6]
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
    set(cb,'YTick',[cy(1),cy(2),cy(end)])
end
set(gcf, 'InvertHardcopy', 'off')
%Save
axis off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS7/' genelist{list} '_CS7_bBG_OFF.png'],'-r800')
h0 = text(5350,1000,{'mRNA level'},'fontweight','FontSize',18,'Color',[1, 1 ,1])
set(h0,'Rotation',90);
text(2900,-350,{'Tb'},'FontSize',12,'Color',[1, 1 ,1])
hold on
plot(linspace(3147,3000,100),linspace(-143,-300,100),'w-','LineWidth',1)
text(5300,600,{'SYS'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(5152,5300,100),linspace(557,600,100),'w-','LineWidth',1)
text(4200,1500,{'EmDisc'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4234,4300,100),linspace(1147,1450,100),'w-','LineWidth',1)
text(2370,480,{'Am'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(2690,2500,100),linspace(577,510,100),'w-','LineWidth',1)
plot(linspace(4843,5400,100),linspace(104,200,100),'w-','LineWidth',1)
text(5400,200,{'ExMes'},'FontSize',12,'Color',[1, 1 ,1])
text(5360,400,{'Stalk'},'FontSize',12,'Color',[1, 1 ,1])
plot(linspace(4673,5304,100),linspace(160,390,100),'w-','LineWidth',1)
text(3700,2430,genelist{list},'FontSize',22,'Color',[1, 1 ,1])
text(3400,2220,{'Carnegie Stage 7'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS7/' genelist{list} '_CS7_bBG_ON.png'],'-r800')

close all
%Pre: CS3
h = figure('visible', 'off');
f1 = patch('Faces',O4.objects(4).data.vertices,'Vertices',O4.vertices,'FaceVertexCData',M7*ones(length(O4.vertices),1),'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',O4.objects(12).data.vertices,'Vertices',O4.vertices,'FaceVertexCData',M8*ones(length(O4.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O4.objects(8).data.vertices,'Vertices',O4.vertices,'FaceVertexCData',M6*ones(length(O4.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',M7*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',M7*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',M5*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
xlim([900,2300])
ylim([700,2100])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gca,'clim',[minV maxV] )
cb=colorbar;
cb.Color = 'w';
cb.FontSize = 16;
cb.Position = cb.Position + [-0.05 0.45 -0 -0.6]
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
    set(cb,'YTick',[cy(1),cy(2),cy(end)])
end
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS3/' genelist{list} '_CS3_bBG_OFF.png'],'-r800')
hold on

h0 = text(2030,1500,{'mRNA level'},'fontweight','FontSize',18,'Color',[1, 1 ,1])
set(h0,'Rotation',90);
text(1300+ DeltaX,1300+ DeltaY,{'ICM'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1326,1325,100)+ DeltaX,linspace(1077,1230,100)+ DeltaY,'w-','LineWidth',1)
text(1500+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1320,1480,100)+ DeltaX,linspace(915,850,100)+ DeltaY,'w-','LineWidth',1)
text(2000+ DeltaX,800+ DeltaY,{'Tb'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1905,2000,100)+ DeltaX,linspace(882,850,100)+ DeltaY,'w-','LineWidth',1)
text(1700+ DeltaX,1300+ DeltaY,{'Epiblast'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1806,1790,100)+ DeltaX,linspace(1124,1240,100)+ DeltaY,'w-','LineWidth',1)
text(2100+ DeltaX,1050+ DeltaY,{'Hypoblast'},'FontSize',22,'Color',[1, 1 ,1])
plot(linspace(1932,2090,100)+ DeltaX,linspace(1067,1065,100)+ DeltaY,'w-','LineWidth',1)
text(1400,2050,genelist{list},'FontSize',22,'Color',[1, 1 ,1])
text(1300,1960,{'Carnegie Stage 3'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS3/' genelist{list} '_CS3_bBG_ON.png'],'-r800')

close all

%Pre CS1-2
h = figure('visible', 'off');
f1 = patch('Faces',O6.objects(4).data.vertices,'Vertices',O6.vertices,'FaceVertexCData',M4*ones(length(O6.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O7.objects(10).data.vertices,'Vertices',O7.vertices,'FaceVertexCData',M3*ones(length(O7.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O8.objects(9).data.vertices,'Vertices',O8.vertices,'FaceVertexCData',M2*ones(length(O8.vertices),1),'FaceColor','interp','LineStyle','none');
f1 = patch('Faces',O9.objects(4).data.vertices,'Vertices',O9.vertices,'FaceVertexCData',M1*ones(length(O9.vertices),1),'FaceColor','interp','LineStyle','none');
xlim([-600,1200])
ylim([200,2000])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gca,'clim',[minV maxV] )
cb=colorbar;
cb.Color = 'w';
cb.FontSize = 16;
cb.Position = cb.Position + [-0.01 0.45 -0 -0.6]
cy=get(cb,'YTick');
try
set(cb,'YTick',[cy(1),cy(3),cy(end)])
catch
    set(cb,'YTick',[cy(1),cy(2),cy(end)])
end
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print('-dpng',['./CS1-2/' genelist{list} '_CS1-2_bBG_OFF.png'],'-r800')
hold on

h0 = text(990,1230,{'mRNA level'},'fontweight','FontSize',18,'Color',[1, 1 ,1])
set(h0,'Rotation',90);

text(-460+ DeltaX,1300+ DeltaY,{'Zygote'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(-310,-360,100)+ DeltaX,linspace(1016,1230,100)+ DeltaY,'w-','LineWidth',1)
text(-20+ DeltaX,1300+ DeltaY,{'4-cell'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(65,15,100)+ DeltaX,linspace(1030,1230,100)+ DeltaY,'w-','LineWidth',1)
text(350+ DeltaX,1300+ DeltaY,{'8-cell'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(472,500,100)+ DeltaX,linspace(1065,1230,100)+ DeltaY,'w-','LineWidth',1)
text(700+ DeltaX,1300+ DeltaY,{'Morula'},'FontSize',22,'Color',[1, 1 ,1])
%plot(linspace(846,906,100)+ DeltaX,linspace(1095,1220,100)+ DeltaY,'w-','LineWidth',1)
text(120,1950,genelist{list},'FontSize',22,'Color',[1, 1 ,1])
text(-50,1850,{'Carnegie Stages 1-2'},'fontweight','bold','FontSize',22,'Color',[1, 1 ,1])
print('-dpng',['./CS1-2/' genelist{list} '_CS1-2_bBG_ON.png'],'-r800')

catch
end

end
