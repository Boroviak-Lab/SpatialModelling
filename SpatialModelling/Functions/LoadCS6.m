function [OBJ2A,sections] = LoadCS6(sections)

%237,243,246,252,,
if strcmp(sections,'3D')==1

load('../Data/SpatialData/CS6_w_CutExMes4.mat')
%1-exmes, 7
%7-exmes bottom, 12
%13 - emdisc,16
%17-pgc, 20
%21-sys, 24
%25,ve,28
%29-am, 32
%33 stqlk 36
%37-tb, 40


%exmes = load('./Data/CutExMes_3.mat')

%exe = load('./Data/CS6_embryomodel_260421_HIGHRES.mat')
%1,emdisc,4
%5,pgc,8
%9 sys,12
%13,ve 16
%17am, 20
%21 stalk, 24
%25 tb, 28
%29 exmes
%32

%OBJ2A.objects([1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,31]) = OBJ2A.objects([13,16,17,20,21,24,25,28,29,32,33,36,37,40,1,6]);

temp.objects = OBJ2A.objects([7,12,1,6]);
OBJ2A.objects([1,4,5,8,9,12,13,16,17,20,21,24,25,28,29]) = OBJ2A.objects([13,16,17,20,21,24,25,28,29,32,33,36,37,40,7]);
OBJ2A.objects(32) = temp.objects(2);
OBJ2A.objects(32).data.vertices = [temp.objects(2).data.vertices;temp.objects(4).data.vertices];
%objects(32).data.vertices = [exmes.OBJ2A_alt.objects(6).data.vertices;exmes.OBJ2A_alt.objects(12).data.vertices]

%29,32,33,36
%13 - emdisc,16
%17-pgc, 20
%21-sys, 24
%25,ve,28
%29-am, 32
%33 stqlk 36
%37-tb, 40
%1-exmes, 6
%7-exmes bottom, 12

%1/6
%7/
%n = size(OBJ2A.vertices,1);
%OBJ2A(1).objects(29) = exmes.OBJ2A_alt(1).objects(1);
%OBJ2A(1).objects(32).data.vertices = [exmes.OBJ2A_alt.objects(6).data.vertices;exmes.OBJ2A_alt.objects(12).data.vertices] + n;
%OBJ2A.vertices = [OBJ2A.vertices; exmes.OBJ2A_alt.vertices];

elseif strcmp(sections,'237')==1 %Partial 
%load('./Data/CS6_embryomodel_260421_HIGHRES.mat')
load('../Data/SpatialData/Cross sections/CS6/CS6_revisions_CROSS_237.mat')
OBJ2A = OBJ2A_1;
OBJ2A.objects([1,13,17,21,25,29]) = OBJ2A_1.objects([1,5,9,13,17,21]);
OBJ2A.objects([4,16,20,24,28,32]) = OBJ2A_1.objects([4,8,12,16,20,24]); %Issue 'cos we missing 8 and 12
 
OBJ2A.objects(12).data.vertices = [];
OBJ2A.objects(12).data.normal = []
OBJ2A.objects(8).data.vertices = [];
OBJ2A.objects(8).data.normal = [];
elseif strcmp(sections,'243')==1
load('../Data/SpatialData/Cross sections/CS6/CS6_revisions_CROSS_243.mat')
OBJ2A = OBJ2A_1;
OBJ2A.objects([1,5,9,13,17,21,25,29]) = OBJ2A_1.objects([1,5,17,11,15,19,27,31]);
OBJ2A.objects([4,8,12,16,20,24,28,32]) = OBJ2A_1.objects([4,10,22,14,18,22,30,34]);
%No PGC or SYS
elseif strcmp(sections,'246')==1
    
%keyboard
load('../Data/SpatialData/Cross sections/CS6/CS6_revisions_CROSS_246_PGCfix.mat')
%load('./Data/Cross sections/CS6/CS6_revisions_CROSS_246.mat')
OBJ2A = OBJ2A_1;
%1=emd,4
%5=PGC,8
%9=VE,12
%13,am, 16
%17,stalk,20
%21 sys,24
%25,Tb,28
%29 exmes, 32

OBJ2A.objects([17,20,1,4,5,8,21,24,13,16,9,12,25,28,29,32]) = OBJ2A.objects([13,16,1,4,5,8,17,20,9,12,21,24,25,28,29,32]);
% 
% if sum(strcmp(tissues,'Am'))~=0
% f1 = patch('Faces',OBJ2A.objects(20).data.vertices,'Vertices',  OBJ2A.vertices,'FaceVertexCData',Output.m1,'FaceColor','interp','LineStyle','none');
% if sum(strcmp(tissues,'EmDisc'))~=0
% f2 = patch('Faces',OBJ2A.objects(4).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m2,'FaceColor','interp','LineStyle','none');
% if sum(strcmp(tissues,'PGC'))~=0
% f3 = patch('Faces',OBJ2A.objects(8).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m3,'FaceColor','interp','LineStyle','none');
% if sum(strcmp(tissues,'Stalk'))~=0
% f4 = patch('Faces',OBJ2A.objects(24).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m8,'FaceColor','interp','LineStyle','none');
% if sum(strcmp(tissues,'VE'))~=0
% f5 = patch('Faces',OBJ2A.objects(16).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m4,'FaceColor','interp','LineStyle','none');
% if sum(strcmp(tissues,'SYS'))~=0
% f6 = patch('Faces',OBJ2A.objects(12).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m5,'FaceColor','interp','LineStyle','none');
% if sum(strcmp(tissues,'Tb'))~=0
% f7 = patch('Faces',OBJ2A.objects(28).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m7,'FaceColor','interp','LineStyle','none');
% if sum(strcmp(tissues,'ExMes'))~=0
% f8 = patch('Faces',OBJ2A.objects(32).data.vertices,'Vertices',OBJ2A.vertices,'FaceVertexCData',Output.m6,'FaceColor','interp','LineStyle','none');



%OBJ2A.objects([1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32]) = 
%OBJ2A.objects();

%Am,EmD,PGC,Stalk,VE,SYS,Tb,ExMes
%17,1,5,21,13,9,25,29,20,4,8,24,16,12,28,32
%20,4,8,24,16,12,28,32
%
%15,1,5,19,11,23,27,31,18,4,10,22,14,26,30,34
%
%Order is: 

%OBJ2A.objects([17,1,5,21,13,9,25,29,20,4,8,24,16,12,28,32]) = OBJ2A_1.objects([15,1,5,19,11,23,27,31,18,4,10,22,14,26,30,34]);
%OBJ2A.objects([8,12,16,20,24,28,32]) = OBJ2A_1.objects([10,22,14,18,22,30,34]);
elseif strcmp(sections,'252')==1
load('../Data/SpatialData/Cross sections/CS6/CS6_revisions_CROSS_252.mat')
OBJ2A = OBJ2A_1;
%No PGC
OBJ2A.objects(5)  = [];
OBJ2A.objects(8)  = [];
OBJ2A.objects([1,9,13,17,21,25,29]) = OBJ2A_1.objects([1,17,5,9,13,21,25]);
OBJ2A.objects([4,12,16,20,24,28,32]) = OBJ2A_1.objects([4,20,8,12,16,24,28]);
else
load('../Data/SpatialData/CS6_embryomodel_260421_HIGHRES.mat')
end
%EMD,1 ; PGC, 5; SYS 9; VE 13; Am 17; Stalk 21; Tb 25; Exmes 29
%EmD  1 ; VE 5 . ;Am  9; Stalk 13; Tb 17; ExMes 21  d
%EmD 1, PGC; VE 11 ; Am 15; Stalk 19; SYS 23; Tb 27; ExMes 31 s
%EmD; PGC 5; VE 11; Am 15; Stalk 19; SYS 23; Tb 27; ExMes 31s
%EmD 1; VE 5;  Am 9; Stalk 13; SYS 17; Tb 21; Exmes 25
