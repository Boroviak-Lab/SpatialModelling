addpath(genpath('./Functions'))

[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
[Output] = loadCS5Scaffold(D,Locations,Shots);

%%%%%%Live editing
%SOX2 high, 
R1 = [-121.9,-32.89,71.36]; 
%MIXL high 
R2 = [109.5,30.07,-130.4];
%Lefty high 
R3 = [-88.17,-47.81,197.9];
%Lefty low 
R4 = [128.6, 24.3, -177];
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
hold on
scatter3(R1(1),R1(2),R1(3),125,'ro','filled')
scatter3(R2(1),R2(2),R2(3),125,'go','filled')
scatter3(R3(1),R3(2),R3(3),125,'bo','filled')
scatter3(R4(1),R4(2),R4(3),125,'yo','filled')

%Fir generate the arc for EmDisc AP via middle point 
R1 = [-121.9,-32.89,71.36]; 
R2 = [109.5,30.07,-130.4];
R12 = [27.2600  -52.2000  -29.5200]
interP = cscvn([R1;R12;R2]');


[ x_ad,y_ad,z_ad,kpi ] = keszito( interP );
Line.vertices = [x_ad,y_ad,z_ad];
save('Data/CS5_EmDisc.mat','Line')
%save('Data/CS5_EmDisc.mat','interP')


f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
hold on
plot3(x_ad,y_ad,z_ad,'k-')
scatter3(R1(1),R1(2),R1(3),125,'ro','filled')
scatter3(R2(1),R2(2),R2(3),125,'go','filled')
scatter3(R12(1),R12(2),R12(3),125,'bo','filled')

%Next we work on the VE
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
hold on
scatter3(R3(1),R3(2),R3(3),125,'ro','filled')
scatter3(R4(1),R4(2),R4(3),125,'go','filled')

%R34 = [(R3(1)+R4(1))/2, Y, (R3(3)+R4(3))/2];
%R4 = [128.6, 24.3, -177];

R34_1 = [-30,-145,78];
R34_2 = [45,-149,-46];
R34_3 = [81,-85,-100];

interP = cscvn([R3;R34_1;R34_2;R34_3;R4]');
[ x_ad,y_ad,z_ad,kpi ] = keszito( interP );

f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
hold on
plot3(x_ad,y_ad,z_ad,'k-')

Line.vertices = [x_ad,y_ad,z_ad];
save('Data/CS5_VE.mat','Line')

f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);

%Amnion 
R6_1 = [-124.6,-30.5,70.07];
R6_2 = [-6.848,10,-2.502];
R6_3 = [-78.88,-3,37.57];

%inner: [-122.5,-31.44,66.92]
%[-70.03,-6.418,37.57]
%[3.93,6.184,-2.586]

interP = cscvn([R6_1;R6_3;R6_2]');
[ x_ad,y_ad,z_ad,kpi ] = keszito( interP );
hold on
plot3(x_ad,y_ad,z_ad,'r-')

Line.vertices = [x_ad,y_ad,z_ad];
save('Data/CS5_Am.mat','Line')

figure(1)
%Process the shots onto scaffold
[Output] = MarmosetGP_CS5(D,Output,'SOX2');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,1);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title('SOX2')
colorbar
[Output] = MarmosetGP_CS5(D,Output,'TFAP2A');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,2);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title('TFAP2A')
colorbar
[Output] = MarmosetGP_CS5(D,Output,'TFAP2C');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,3);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title('TFAP2C')
colorbar
[Output] = MarmosetGP_CS5(D,Output,'MIXL1');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,4);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title('MIXL1')
colorbar



figure(2)
%Process the shots onto scaffold
[Output] = MarmosetGP_CS5(D,Output,'SOX2');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,1);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m1),[1 N] );        
r = normrnd(Output.m1(ind),Output.s1(ind));
scatter( (ind+100).*cos(u),(ind+100).*sin(u),20,r,'filled')
title('SOX2')
colorbar
[Output] = MarmosetGP_CS5(D,Output,'TFAP2A');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,2);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m1),[1 N] );        
r = normrnd(Output.m1(ind),Output.s1(ind));
scatter((ind+100).*cos(u),(ind+100).*sin(u),20,r,'filled')
title('TFAP2A')
colorbar
[Output] = MarmosetGP_CS5(D,Output,'TFAP2C');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,3);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m1),[1 N] );        
r = normrnd(Output.m1(ind),Output.s1(ind));
scatter((ind+100).*cos(u),(ind+100).*sin(u),20,r,'filled')
title('TFAP2C')
colorbar
[Output] = MarmosetGP_CS5(D,Output,'MIXL1');
Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS5(Output,Line,'Line');
subplot(2,2,4);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m1),[1 N] );        
r = normrnd(Output.m1(ind),Output.s1(ind));
scatter((ind+100).*cos(u),(ind+100).*sin(u),20,r,'filled')
title('MIXL1')
colorbar




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Next we do CS6 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[OBJ2A,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS6');
[Output] = loadCS6Scaffold(D,Locations,Shots);
PlotEmbryoCS6(OBJ2A,'all');


[Output] = MarmosetGP_CS6(D,Output,'MIXL1');
[Output] = MarmosetGPInfer_CS6(Output,OBJ2A);
PlotEmbryoCS6GP(Output,OBJ2A,{'EmDisc','Stalk'},1);

[Output] = MarmosetGP_CS6(D,Output,'SOX2');
[Output] = MarmosetGPInfer_CS6(Output,OBJ2A);
PlotEmbryoCS6GP(Output,OBJ2A,{'EmDisc','Stalk'},2);


[Output] = MarmosetGP_CS6(D,Output,'LEFTY2');
[Output] = MarmosetGPInfer_CS6(Output,OBJ2A);
PlotEmbryoCS6GP(Output,OBJ2A,{'EmDisc','VE','Stalk'},3);
PlotEmbryoCS6GP(Output,OBJ2A,{'EmDisc','VE'},4);

%SOX2 high
R1 = [-160.7,-255.6,29.87];
%Mixl high
R2 = [117.4,261.2,-34.65];
R12_1 = [-65,-50,14]
R12_2 = [112,94.5,6.265],[80.21,101.2,6.116]
R12_3 = [120,176,-8]

interP = cscvn([R1;R12_1;R12_2;R12_3;R2]');
[ x_ad,y_ad,z_ad,kpi ] = keszito( interP );

Line.vertices = [x_ad,y_ad,z_ad];
save('Data/CS6_EmDisc.mat','Line')

%save('Data/CS6_EmDisc.mat','interP')

f2 = patch('Faces',OBJ2A.objects(4).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2A.objects(24).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
hold on
scatter3(R1(1),R1(2),R1(3),125,'ro','filled')
scatter3(R2(1),R2(2),R2(3),125,'go','filled')
plot3(x_ad,y_ad,z_ad,'k-')


%Lefty high
R3 = [-24.84,-127.4,134.5];
%Lefty low
R4 = [207.1,156.2,8.686];
R34_1 = [80,-25,85];
R34_2 = [200,50,40];

interP = cscvn([R3;R34_1;R34_2;R4]');
[ x_ad,y_ad,z_ad,kpi ] = keszito( interP );
plot3(x_ad,y_ad,z_ad,'k-')

Line.vertices = [x_ad,y_ad,z_ad];
save('Data/CS6_VE.mat','Line')

f5 = patch('Faces',OBJ2A.objects(16).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f2 = patch('Faces',OBJ2A.objects(4).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2A.objects(24).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
hold on
scatter3(R3(1),R3(2),R3(3),125,'ro','filled')
scatter3(R4(1),R4(2),R4(3),125,'go','filled')

f5 = patch('Faces',OBJ2A.objects(16).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2A.objects(4).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2A.objects(24).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
hold on
scatter3(R1(1),R1(2),R1(3),125,'ro','filled')
scatter3(R2(1),R2(2),R2(3),125,'go','filled')

f1 = patch('Faces',OBJ2A.objects(20).data.vertices,'Vertices',  OBJ2A.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2A.objects(4).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2A.objects(8).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2A.objects(24).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2A.objects(16).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now do CS7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[OBJ3,section] = LoadCS7('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS7');
[Output] = loadCS7Scaffold(D,Locations,Shots);
%PlotEmbryoCS6(OBJ1,tissues);


[Output] = MarmosetGP_CS7(D,Output,'MIXL1');
[Output] = MarmosetGPInfer_CS7(Output,OBJ3);
PlotEmbryoCS7GP(Output,OBJ3,{'EmDisc','Stalk'},1);

[Output] = MarmosetGP_CS7(D,Output,'SOX2');
[Output] = MarmosetGPInfer_CS7(Output,OBJ3);
PlotEmbryoCS7GP(Output,OBJ3,{'EmDisc','Stalk'},2);

%SOX2 high 
R1 = [-1163,-537,27.09,]
%Mixl1 high 
R2 = [23.15,-285.8,33]
R12_1 = [-890,-405,90]
R12_2 = [-665,-150,50]
R12_3 = [-252,-190,-5] 

interP = cscvn([R1;R12_1;R12_2;R12_3;R2]');
[ x_ad,y_ad,z_ad,kpi ] = keszito( interP );
plot3(x_ad,y_ad,z_ad,'k-')

Line.vertices = [x_ad,y_ad,z_ad];
save('Data/CS7_EmDisc.mat','Line')
%save('Data/CS7_EmDisc.mat','interP')

f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[2,51,191]/255,'LineStyle','none','FaceAlpha',.1);
f8 = patch('Faces',OBJ3.objects(30).data.vertices,'Vertices',OBJ3.vertices,'FaceColor',[96, 56, 19]/255,'LineStyle','none','FaceAlpha',.1);
hold on
scatter3(R1(1),R1(2),R1(3),125,'ro','filled')
scatter3(R2(1),R2(2),R2(3),125,'go','filled')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now do CS6_2
[OBJ1,section] = LoadCS6_2('3D');
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS62');
[Output] = loadCS6Scaffold2(D,Locations,Shots);

[Output] = MarmosetGP_CS62(D,Output,'SOX2');
[Output] = MarmosetGPInfer_CS6(Output,OBJ1);
PlotEmbryoCS6GP_2(Output,OBJ1,{'EmDisc','Stalk'},1);



[Output] = MarmosetGP_CS62(D,Output,'MIXL1');
[Output] = MarmosetGPInfer_CS6(Output,OBJ1);
PlotEmbryoCS6GP_2(Output,OBJ1,{'EmDisc','Stalk'},2);



%SOX2 high
R1 = ([40.85,-133.9,86.03] + [-11.95,-154.1,103] )/2
%Mixl high
R2 = ([-45.79,148,-120.2] + [-16.98,-139.3,-120.4]) /2
R12 = ([28.94,-2.277,33.91] + [-15.69,-53.11,66.11])/2
%R1 = [-160.7,-255.6,29.87];
%Mixl high
%R2 = [117.4,261.2,-34.65];
%R12_1 = [-65,-50,14]
%R12_2 = [112,94.5,6.265],[80.21,101.2,6.116]
%R12_3 = [120,176,-8]

interP = cscvn([R1;R12;R2]');
[ x_ad,y_ad,z_ad,kpi ] = keszito( interP );

Line.vertices = [x_ad,y_ad,z_ad];
save('Data/CS62_EmDisc.mat','Line')
PlotEmbryoCS6(Output,OBJ1,{'EmDisc','Stalk'},2);


