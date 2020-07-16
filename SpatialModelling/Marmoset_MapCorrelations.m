addpath(genpath('./Functions/'))
D3 = importdata('Data/Keycorrect.csv')
Dexp = importdata('Data/NormData.csv');

[DexpCS5,Arc,Xtest,OutputCS5] = loadCS5Probs(D3,Dexp);
[DexpCS6,Arc,Xtest,OutputCS6] = loadCS6Probs(D3,Dexp);

[Output1] = Marmoset3D_CS5_wholesurface(DexpCS5,OutputCS5,'Base2','SOX2')
[Output2] = Marmoset3D_CS6_wholesurface(DexpCS6,OutputCS6,'Base2','SOX2')
%[Output3] = Marmoset3D_CS7_wholesurface(DexpCS7,OutputCS7,'Base2','SOX2')
max1 = max([Output1.m_1;Output2.m_1;Output3.m_1]);
min1 = min([Output1.m_1;Output2.m_1;Output3.m_1]);

[Output4] = Marmoset3D_CS5_wholesurface(DexpCS5,OutputCS5,'Base2','HESX1')
[Output5] = Marmoset3D_CS6_wholesurface(DexpCS6,OutputCS6,'Base2','HESX1')
%[Output6] = Marmoset3D_CS7_wholesurface(DexpCS7,OutputCS7,'Base2','HESX1')
max2 = max([Output4.m_1;Output4.m_1;Output6.m_1]);
min2 = min([Output4.m_1;Output4.m_1;Output6.m_1]);

% Sample values 
h = 10; % height 
ra = 20; % radius
tht = linspace(0,2*pi,100); z = linspace(0,h,20);
xa = repmat(ra*cos(tht),20,1); ya = repmat(ra*sin(tht),20,1); 
za = repmat(z',1,100);

% To close the ends 
X = [xa*0; flipud(xa); (xa(1,:))*0]; Y = [ya*0; flipud(ya); (ya(1,:))*0]; 
Z = [za; flipud(za); za(1,:)];

% Draw cylinder 
[TRI,v]= surf2patch(X,Y,Z,'triangle'); 
V = v;
Y = TRI;
for i = 1:30
    v2 = v;
    v2(:,3) = v2(:,3) + (i*10)*ones(size(v,1),1);
    V = [V;v2];
    m = max(max(Y));
    m1 = TRI+ m;
    Y = [Y;m1];
end

patch('Vertices',V,'Faces',Y,'FaceVertexCData',C,'FaceColor','interp','LineStyle','none')

C1 = linspace(0.5,0.7,size(V,1))';
patch('Vertices',V,'Faces',Y,'FaceVertexCData',C1,'FaceColor','interp','LineStyle','none')
%V1 = V; V1(:,[1,3]) = V1(:,[3,1]); 
V1 = V; V1(:,[2,3]) = V1(:,[3,2]); 
V1(:,2) = V1(:,2) -200;
V1(:,1) = V1(:,1) +200;

%[Output2] = Marmoset3D_CS5_wholesurface(DexpCS5,OutputCS5,'Base2','POU5F1')
%[Output3] = Marmoset3D_CS5_wholesurface(DexpCS5,OutputCS5,'Base2','T')
%[Output4] = Marmoset3D_CS5_wholesurface(DexpCS5,OutputCS5,'Base2','MIXL1')

%OBJ1.vertices = Quaternion3(rotx,[0,0,1],OBJ1.vertices); 
x = [0 0 0;
    1 0 0;
    1 1 0;
    0 1 0;
    0 0 1;
    1 0 1;
    1 1 1;
    0 1 1];
y = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

%x = [1 1 2 2];
%y = [0 1 1 0];

%x =[0, 0 ; 0, 1; 1, 1; 1,0];
%y = y;%[1,2,3,4];
X = x;
Y = y;
for i = 1:30    
X = [X; x + [i,0,0;i,0,0;i,0,0;i,0,0;i,0,0;i,0,0;i,0,0;i,0,0] ];       
%X = [x+i*[1 1 1 1]];
Y = [Y;y+(i)*8];
end
C = linspace(0.5,0.7,size(X,1))';
patch('Faces',Y,'Vertices',X,'FaceVertexCData',C,'FaceColor','interp','LineStyle','none')
camlight('headlight')

Xt = X; Xt(:,1) = Xt(:,1)*6 -100; Xt(:,2) = Xt(:,2)*50 + 100; Xt(:,3)=Xt(:,3)*6;
Xt(:,2) = Xt(:,2)+50;
Xt(:,[1:2]) = Xt(:,[2,1]);

%x1=x; x1(:,1) = x1(:,1)+1;
%patch('Vertices',x1,'Faces',y,...
 %     'FaceVertexCData',hsv(6),'FaceColor','flat')
%patch(x,y,'red')
bronze=[178/255,86/255,13/255];
silver=[197/255,206/255,212/255];
%gold=[249/255,166/255,2/255];

gold=[215/255,190/255,105/255];
jade=[124/255,157/255,142/255];
ruby=[224/255,17/255,95/255];


n = 100;                %// number of colors
cmap(1,:) = [0,0,0];   %// color first row - red
cmap(2,:) = bronze;   %// color 25th row - green
cmap(3,:) = [232,224,9]/255;   %// color 50th row - blue
[X1,Y1] = meshgrid([1:3],[1:100]);  %// mesh of indices
chriscmap = interp2(X1([1,25,100],:),Y1([1,25,100],:),cmap,X1,Y1); %// interpolate colormap
colormap(chriscmap) %// set color map

%Amn, Em, 
OBJ1=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS5.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ1.vertices = Quaternion3(20,[1,0,0],OBJ1.vertices); 

h = figure(1)
subplot(2,3,1);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{1,1}+OutputCS5.Convf{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,1}+OutputCS5.Convf{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{3,1}+OutputCS5.Convf{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{4,1}+OutputCS5.Convf{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{5,1}+OutputCS5.Convf{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{6,1}+OutputCS5.Convf{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')

subplot(2,3,2);
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,1},'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',Y,'Vertices',Xt,'FaceVertexCData',C,'FaceColor','interp','LineStyle','none')
f3 = patch('Vertices',V1,'Faces',Y,'FaceVertexCData',C1,'FaceColor','interp','LineStyle','none')
view([74.2257,61.0724])
axis equal
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight('headlight')

subplot(2,3,4);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')

subplot(2,3,5);
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,1},'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',Y,'Vertices',Xt,'FaceVertexCData',C,'FaceColor','interp','LineStyle','none')
%f3 = patch('Vertices',V,'Faces',Y,'FaceVertexCData',C1,'FaceColor','interp','LineStyle','none')
f3 = patch('Vertices',V1,'Faces',Y,'FaceVertexCData',C1,'FaceColor','interp','LineStyle','none')
view([74.2257,61.0724])
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight('headlight')

subplot(2,3,3);
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_1,'FaceColor','interp','LineStyle','none');
view([74.2257,61.0724])
%hold on
title('SOX2')
set(gca,'clim',[min1 max1] )
axis equal
axis off
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)

%x = [0 0 1 1];
%y = [0 1 1 0];


%camlight
subplot(2,3,6);
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output4.m_1,'FaceColor','interp','LineStyle','none');

%f3 = patch('Faces',Y,'Vertices',X*200,'FaceVertexCData',C,'FaceColor','interp','LineStyle','none')

view([74.2257,61.0724])
%hold on
title('HESX1')
set(gca,'clim',[min2 max2] )
axis equal
axis off
%camlight
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)

print(h,['~/Desktop/ESCs_CS5_sox2.png'],'-dpng','-r2100'); 
clf


%%%%%

h = figure(1)
subplot(2,3,1);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{1,1}-OutputCS5.Convf{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,1}-OutputCS5.Convf{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{3,1}-OutputCS5.Convf{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{4,1}-OutputCS5.Convf{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{5,1}-OutputCS5.Convf{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{6,1}-OutputCS5.Convf{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')

subplot(2,3,2);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')


subplot(2,3,3);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{1,1}+OutputCS5.Convf{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,1}+OutputCS5.Convf{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{3,1}+OutputCS5.Convf{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{4,1}+OutputCS5.Convf{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{5,1}+OutputCS5.Convf{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{6,1}+OutputCS5.Convf{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')


subplot(2,3,4);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{1,1}-OutputCS5.Naivef{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,1}-OutputCS5.Naivef{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{3,1}-OutputCS5.Naivef{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{4,1}-OutputCS5.Naivef{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{5,1}-OutputCS5.Naivef{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{6,1}-OutputCS5.Naivef{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')


subplot(2,3,5);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')


subplot(2,3,6);
f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{1,1}+OutputCS5.Naivef{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,1}+OutputCS5.Naivef{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{3,1}+OutputCS5.Naivef{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{4,1}+OutputCS5.Naivef{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{5,1}+OutputCS5.Naivef{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{6,1}+OutputCS5.Naivef{6,1},'FaceColor','interp','LineStyle','none');
view([-45.9432,39.0795])
axis equal
axis off
set(gca,'clim',[.5 .7] )
camlight('headlight')


h = figure(2)
% subplot(2,3,1);
% f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{1,2},'FaceColor','interp','LineStyle','none');
% hold on
% f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,2},'FaceColor','interp','LineStyle','none');
% f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{3,2},'FaceColor','interp','LineStyle','none');
% f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{4,2},'FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{5,2},'FaceColor','interp','LineStyle','none');
% f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{6,2},'FaceColor','interp','LineStyle','none');
% view([-45.9432,39.0795])
% axis equal
% axis off
% %set(gca,'clim',[.2 .8] )
% camlight('headlight')

% subplot(2,3,2);
% f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,3},'FaceColor','interp','LineStyle','none');
% view([74.2257,61.0724])
% axis equal
% set(gca,'clim',[.78 .82] )
% axis equal
% axis off
% camlight('headlight')
% cb=colorbar;
% cb.Position = cb.Position + [0.1 -0.02 0 0]
% cy=get(cb,'YTick');
% set(cb,'YTick',[cy(1),cy(3),cy(end)])
% set(gca,'fontsize', 12)

subplot(2,3,2);
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Conv{2,4},'FaceColor','interp','LineStyle','none');
view([74.2257,61.0724])
title('SOX2')
axis equal
set(gca,'clim',[.76 .87] )
axis equal
axis off
camlight('headlight')
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(5),cy(end)])
set(gca,'fontsize', 12)

% subplot(2,3,4);
% f1 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{1,2},'FaceColor','interp','LineStyle','none');
% hold on
% f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,2},'FaceColor','interp','LineStyle','none');
% f3 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{3,2},'FaceColor','interp','LineStyle','none');
% f4 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{4,2},'FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{5,2},'FaceColor','interp','LineStyle','none');
% f6 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{6,2},'FaceColor','interp','LineStyle','none');
% view([-45.9432,39.0795])
% axis equal
% axis off
% %set(gca,'clim',[.2 .8] )
% camlight('headlight')

% subplot(2,3,5);
% f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,3},'FaceColor','interp','LineStyle','none');
% view([74.2257,61.0724])
% set(gca,'clim',[.78 .82] )
% axis equal
% axis off
% camlight('headlight')
% cb=colorbar;
% cb.Position = cb.Position + [0.1 -0.02 0 0]
% cy=get(cb,'YTick');
% set(cb,'YTick',[cy(1),cy(3),cy(end)])
% set(gca,'fontsize', 12)

subplot(2,3,5);
f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive{2,4},'FaceColor','interp','LineStyle','none');
view([74.2257,61.0724])

set(gca,'clim',[.76 .87] )
axis equal
axis off
camlight('headlight')
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(5),cy(end)])
set(gca,'fontsize', 12)

% subplot(2,3,3);
% f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output1.m_1,'FaceColor','interp','LineStyle','none');
% view([74.2257,61.0724])
% %hold on
% set(gca,'clim',[min1 max1] )
% axis equal
% axis off
% %camlight
% subplot(2,3,6);
% f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output4.m_1,'FaceColor','interp','LineStyle','none');
% view([74.2257,61.0724])
% %hold on
% set(gca,'clim',[min2 max2] )
% axis equal
% axis off
% %camlight

h = figure(3)
print(h,['~/Desktop/ESCs_CS5_mean.png'],'-dpng','-r2100'); 
clf

%subplot(1,3,2);
%f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',OutputCS5.Naive,'FaceColor','interp','LineStyle','none');
%view([9 15])
%set(gca,'clim',[0.4 0.7] )
%hold on
%f2 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output.m_0,'FaceColor','interp','LineStyle','none');
%camlight
%axis off

%Xex2 = [OBJ2.vertices(OBJ2.objects(28).data.vertices,:)]; %OBJ2.vertices(OBJ2.objects(32).data.vertices,:)
%Xem2 = [OBJ2.vertices(OBJ2.objects(8).data.vertices,:)];
%Xtroph2 = [OBJ2.vertices(OBJ2.objects(20).data.vertices,:)];
%Xam2 = [OBJ2.vertices(OBJ2.objects(4).data.vertices,:)];
%Xsys2 = [OBJ2.vertices(OBJ2.objects(16).data.vertices,:)];
%Xve2 = [OBJ2.vertices(OBJ2.objects(24).data.vertices,:)];
%Xpgc2 = [OBJ2.vertices(OBJ2.objects(12).data.vertices,:)];

OBJ2=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS6.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
%OBJ2.vertices = Quaternion3(10,[0,0,1],OBJ2.vertices); 
OBJ2.vertices = Quaternion3(10,[0,1,0],OBJ2.vertices); 
%OBJ2.vertices = Quaternion3(10,[0,1,0],OBJ2.vertices);
h = figure(2)
subplot(2,3,1);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{6,1},'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{7,1},'FaceColor','interp','LineStyle','none');
view([-216.0169,-87.7009])
%hold on
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight('headlight')



subplot(2,3,2);
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{2},'FaceColor','interp','LineStyle','none');
view([-202.3130,-70.8686])
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight('headlight')

%view([-208.5981,-51.2960])
%view([-202.3130,-70.8686])

subplot(2,3,4);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{1,1},'FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{2,1},'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{3,1},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{4,1},'FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{5,1},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{6,1},'FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{7,1},'FaceColor','interp','LineStyle','none');
view([-216.0169,-87.7009])
%hold on
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight('headlight')


subplot(2,3,5);
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{2,1},'FaceColor','interp','LineStyle','none');
view([-202.3130,-70.8686])
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight('headlight')



subplot(2,3,3);
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',Output2.m_1,'FaceColor','interp','LineStyle','none');
view([-202.3130,-70.8686])
%hold on
title('SOX2')
set(gca,'clim',[min1 max1] )
axis equal
axis off
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)
%camlight


subplot(2,3,6);
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',Output5.m_1,'FaceColor','interp','LineStyle','none');
view([-202.3130,-70.8686])
title('HESX1')
%hold on
set(gca,'clim',[min2 max2] )
axis equal
axis off
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)
%camlight

print(h,['~/Desktop/ESCs_CS6_sox2.png'],'-dpng','-r2100'); 





h = figure(3)

% subplot(2,3,2);
% f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{2,2},'FaceColor','interp','LineStyle','none');
% view([-202.3130,-70.8686])
% %set(gca,'clim',[.2 .8] )
% axis equal
% axis off
% camlight('right')
% set(gca,'clim',[.78 .82] )
% axis equal
% axis off
% camlight('headlight')
% cb=colorbar;
% cb.Position = cb.Position + [0.1 -0.02 0 0]
% cy=get(cb,'YTick');
% set(cb,'YTick',[cy(1),cy(3),cy(end)])
% set(gca,'fontsize', 12)

subplot(2,3,3);
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Conv{2,4},'FaceColor','interp','LineStyle','none');
view([-202.3130,-70.8686])
%set(gca,'clim',[.2 .8] )
axis equal
axis off
camlight('right')
set(gca,'clim',[.76 .87] )
axis equal
axis off
camlight('headlight')
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(5),cy(end)])
set(gca,'fontsize', 12)

% subplot(2,3,5);
% f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{2,3},'FaceColor','interp','LineStyle','none');
% view([-202.3130,-70.8686])
% %set(gca,'clim',[.2 .8] )
% axis equal
% axis off
% set(gca,'clim',[.78 .82] )
% axis equal
% axis off
% camlight('headlight')
% cb=colorbar;
% cb.Position = cb.Position + [0.1 -0.02 0 0]
% cy=get(cb,'YTick');
% set(cb,'YTick',[cy(1),cy(3),cy(end)])
% set(gca,'fontsize', 12)

subplot(2,3,6);
f2 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',OutputCS6.Naive{2,4},'FaceColor','interp','LineStyle','none');
view([-202.3130,-70.8686])
%set(gca,'clim',[.2 .8] )
axis equal
axis off
camlight('headlight')
set(gca,'clim',[.76 .87] )
axis equal
axis off
camlight('headlight')
cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(5),cy(end)])
set(gca,'fontsize', 12)




print(h,['~/Desktop/ESCs_CS6_mean_oct4.png'],'-dpng','-r2100'); 


%CS5  view([-45.9432,39.0795])
%view([74.2257,61.0724])
%alt view -58.4543,15.5388


%Could we try CS7?
OBJ3=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS7.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ4=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS7_section_cut.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
%OBJ2.vertices = Quaternion3(10,[0,0,1],OBJ2.vertices); 
%OBJ2.vertices = Quaternion3(10,[0,1,0],OBJ2.vertices); 
subplot(2,3,1);
f1 = patch('Faces',OBJ4.objects(6).data.vertices,'Vertices',OBJ4.vertices,'FaceVertexCData',OutputCS7.Conv{1,2},'FaceColor','interp','LineStyle','none');
hold on
f3 = patch('Faces',OBJ4.objects(10).data.vertices,'Vertices',OBJ4.vertices,'FaceVertexCData',OutputCS7.Conv{3,2},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ3.objects(30).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Conv{4,1},'FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ3.objects(24).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS6.Conv{5},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Conv{5,1},'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Conv{2,1},'FaceColor','interp','LineStyle','none');
view([200,-3])
%hold on
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight

subplot(2,3,2);
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Conv{2,1},'FaceColor','interp','LineStyle','none');
view([200,-3])
%hold on
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight



subplot(2,3,3);
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Conv{2,1},'FaceColor','interp','LineStyle','none');
view([200,-3])
%hold on
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight


subplot(2,3,4);
f1 = patch('Faces',OBJ4.objects(6).data.vertices,'Vertices',OBJ4.vertices,'FaceVertexCData',OutputCS7.Naive{1,2},'FaceColor','interp','LineStyle','none');
hold on
f3 = patch('Faces',OBJ4.objects(10).data.vertices,'Vertices',OBJ4.vertices,'FaceVertexCData',OutputCS7.Naive{3,2},'FaceColor','interp','LineStyle','none');
f4 = patch('Faces',OBJ3.objects(30).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Naive{4,1},'FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ3.objects(24).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS6.Conv{5},'FaceColor','interp','LineStyle','none');
f6 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Naive{5,1},'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Naive{2,1},'FaceColor','interp','LineStyle','none');
view([200,-3])
%hold on
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight

subplot(2,3,5);
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',OutputCS7.Naive{2,1},'FaceColor','interp','LineStyle','none');
view([200,-3])
%hold on
set(gca,'clim',[.5 .7] )
axis equal
axis off
camlight


subplot(2,3,3);
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',Output3.m_1,'FaceColor','interp','LineStyle','none');
view([200,-3])
%hold on
set(gca,'clim',[min1 max1] )
axis equal
axis off
camlight

subplot(2,3,6);
f2 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',Output6.m_1,'FaceColor','interp','LineStyle','none');
view([200,-3])
%hold on
set(gca,'clim',[min2 max2] )
axis equal
axis off
camlight



%Now for blast

%Dexp2 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Final_AllGoodShots_wCS6r/NormData2.csv');


Dcorr = importdata("~/Desktop/Thorsten/FINAL/CCAFigForPaper_wCS6_v2//DimRed/ESC_corr_wCS6.csv")
targs = Dcorr.textdata(1,2:end);
targs = strrep(targs,'"','');
targs = strrep(targs,'RNA.','');
escs = Dcorr.textdata(2:end,1);


idnss0a = find(contains(targs,'Tb_CS3'));
idnss0 = find(contains(targs,'Hyp_CS3'));
idnss = find(contains(targs,'Epi_CS3'));
idnss1 = find(contains(escs,'ESC_primed'));
idnss2 = find(contains(escs,'EMS3_PLAXA'));




R1 = Dcorr.data(idnss1,idnss); %primes vs epi
R2 = Dcorr.data(idnss2,idnss); %plaxa vs epi

R3 = Dcorr.data(idnss1,idnss0); %primed vs hyp
R4 = Dcorr.data(idnss2,idnss0); %plaxa vs hyp

R5 = Dcorr.data(idnss1,idnss0a); %primed vs tb
R6 = Dcorr.data(idnss2,idnss0a); %plaxa vs tb

OBJ0a=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/Blastocyst_halfcut_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ0b=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/Blastocyst_slice_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')


%OBJ0=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/blastocyst.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
h = figure(3)
subplot(2,3,1);
f1 = patch('Faces',OBJ0a.objects(12).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R5))./(size(R5,1)*size(R5,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0a.objects(8).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R3))./(size(R3,1)*size(R3,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0a.objects(4).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R1))./(size(R1,1)*size(R1,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(12).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
    view([132.7748,-12.8427])
set(gca,'clim',[.5 .7] )
%cb=colorbar;
%cb.Position = cb.Position + [0.1 0.15 -0 -0.3]
%cy=get(cb,'YTick');
%set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[1,0,0],'LineStyle','none','FaceAlpha',0.1);

subplot(2,3,4);
f1 = patch('Faces',OBJ0a.objects(12).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R6))./(size(R6,1)*size(R6,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0a.objects(8).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R4))./(size(R4,1)*size(R4,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0a.objects(4).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R2))./(size(R2,1)*size(R2,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
    view([132.7748,-12.8427])
    set(gca,'clim',[.5 .7] )
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight



%OBJ0=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/blastocyst.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
subplot(2,3,2);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R5))./(size(R5,1)*size(R5,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R3))./(size(R3,1)*size(R3,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R1))./(size(R1,1)*size(R1,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(12).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([164.2550,-40.2050])
set(gca,'clim',[.5 .7] )
%cb=colorbar;
%cb.Position = cb.Position + [0.1 0.15 -0 -0.3]
%cy=get(cb,'YTick');
%set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[1,0,0],'LineStyle','none','FaceAlpha',0.1);

subplot(2,3,5);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R6))./(size(R6,1)*size(R6,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R4))./(size(R4,1)*size(R4,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R2))./(size(R2,1)*size(R2,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([164.2550,-40.2050])
set(gca,'clim',[.5 .7] )
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight






Dexp2 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Final_AllGoodShots_wCS6r/NormData3.csv');

targs = Dexp2.textdata(1,2:end);
escs = Dexp2.textdata(2:end,1);

idnss0 = find(contains(targs,'Tb_CS3'));
idnss1 = find(contains(targs,'Hyp_CS3'));
idnss2 = find(contains(targs,'Epi_CS3'));

id1 = find(contains(escs,'SOX2'));
id2 = find(contains(escs,'T'));
id3 = find(contains(escs,'POU5F1'));
id4 = find(contains(escs,'HESX1'));

R6 = Dexp2.data(id3,idnss0);
R4 = Dexp2.data(id3,idnss1);
R2 = Dexp2.data(id3,idnss2);

subplot(2,3,3);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R6))./(size(R6,1)*size(R6,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R4))./(size(R4,1)*size(R4,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R2))./(size(R2,1)*size(R2,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([164.2550,-40.2050])
%set(gca,'clim',[.2 .8] )
set(gca,'clim',[min2 max2] )
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
%camlight


R6 = Dexp2.data(id2,idnss0);
R4 = Dexp2.data(id2,idnss1);
R2 = Dexp2.data(id2,idnss2);


subplot(2,3,6);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R6))./(size(R6,1)*size(R6,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R4))./(size(R4,1)*size(R4,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R2))./(size(R2,1)*size(R2,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([164.2550,-40.2050])
%set(gca,'clim',[.2 .8] )
set(gca,'clim',[min2 max2] )
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
%camlight





print(h,['~/Desktop/ESCs_CS3.png'],'-dpng','-r1500'); 




%Dexp 
% indC = find(strcmp(Dexp.textdata(2:end,1),'SOX2'));
% 
% indA = D3.textdata(find(contains(D3.textdata(:,12),'Hyp_CS3')),1);
% indB = D3.textdata(find(contains(D3.textdata(:,12),'Epi_CS3')),1);
% indD = Dexp.textdata(1,2:end);
% for ii = 1:length(indA)
%     find(strcmp(indA{ii},indD)==1)
% end
% 
% %M1 = mean(Dexp.data(indC,indA));
% %M2 = mean(Dexp.data(indC,indB));
% 
% %for i = 1:length(Dexp.textdata(1,2:end)
% 
% %D3 = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/Keycorrect_CPAll2.csv')
% 
% 
% f1 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output.m_0,'FaceColor','interp','LineStyle','none');
% hold on
% f2 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',Output.m_0,'FaceColor','interp','LineStyle','none');
% f3 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',[220,220,220]/255,'LineStyle','none','FaceAlpha',0.1);
% f4 = patch('Faces',OBJ1.objects(28).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',[220,220,220]/255,'LineStyle','none','FaceAlpha',0.1);
% f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',[220,220,220]/255,'LineStyle','none','FaceAlpha',0.1);
% f7 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',[220,220,220]/255,'LineStyle','none','FaceAlpha',0.1);
% f8 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',[220,220,220]/255,'LineStyle','none','FaceAlpha',0.1);
% xlim([-600 600])







Dcorr = importdata("~/Desktop/Thorsten/FINAL/CCAFigForPaper_wCS6_v2//DimRed/ESC_corr_wCS6.csv")
targs = Dcorr.textdata(1,2:end);
targs = strrep(targs,'"','');
targs = strrep(targs,'RNA.','');
escs = Dcorr.textdata(2:end,1);
idnss0a = find(contains(targs,'Tb_CS3'));
idnss0 = find(contains(targs,'Hyp_CS3_2'));
idnss = find(contains(targs,'Epi_CS3_2'));
idnss1 = find(contains(escs,'ESC_primed'));
idnss2 = find(contains(escs,'EMS3_PLAXA'));


R1 = Dcorr.data(idnss1,idnss); %primes vs epi
R2 = Dcorr.data(idnss2,idnss); %plaxa vs epi

R3 = Dcorr.data(idnss1,idnss0); %primed vs hyp
R4 = Dcorr.data(idnss2,idnss0); %plaxa vs hyp

R5 = Dcorr.data(idnss1,idnss0a); %primed vs tb
R6 = Dcorr.data(idnss2,idnss0a); %plaxa vs tb

OBJ0a=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/Blastocyst_halfcut.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ0b=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/Blastocyst_slice.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')


%OBJ0=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/blastocyst.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
subplot(2,3,1);
f1 = patch('Faces',OBJ0a.objects(12).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R5))./(size(R5,1)*size(R5,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0a.objects(8).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R3))./(size(R3,1)*size(R3,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0a.objects(4).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R1))./(size(R1,1)*size(R1,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(12).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([159.6675,9.3268])
set(gca,'clim',[.2 .8] )
%cb=colorbar;
%cb.Position = cb.Position + [0.1 0.15 -0 -0.3]
%cy=get(cb,'YTick');
%set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[1,0,0],'LineStyle','none','FaceAlpha',0.1);

subplot(2,3,4);
f1 = patch('Faces',OBJ0a.objects(12).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R6))./(size(R6,1)*size(R6,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0a.objects(8).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R4))./(size(R4,1)*size(R4,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0a.objects(4).data.vertices,'Vertices',OBJ0a.vertices,'FaceVertexCData',sum(sum(R2))./(size(R2,1)*size(R2,2))*ones(length(OBJ0a.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([132.7748,-12.8427])
set(gca,'clim',[.2 .8] )
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight



%OBJ0=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/blastocyst.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
subplot(2,3,2);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R5))./(size(R5,1)*size(R5,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R3))./(size(R3,1)*size(R3,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R1))./(size(R1,1)*size(R1,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(12).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([164.2550,-40.2050])
set(gca,'clim',[.2 .8] )
%cb=colorbar;
%cb.Position = cb.Position + [0.1 0.15 -0 -0.3]
%cy=get(cb,'YTick');
%set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[1,0,0],'LineStyle','none','FaceAlpha',0.1);

subplot(2,3,5);
f1 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R6))./(size(R6,1)*size(R6,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
hold on
f0 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R4))./(size(R4,1)*size(R4,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',sum(sum(R2))./(size(R2,1)*size(R2,2))*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%patch('Faces',OBJ0.objects(4).data.vertices,'Vertices',OBJ0.vertices,'FaceColor',[.1,.1,.1],'LineStyle','none','FaceAlpha',0.1);
view([164.2550,-40.2050])
set(gca,'clim',[.2 .8] )
set(gca,'fontsize', 12)
set(gcf,'color','w');
axis equal
axis off
camlight

