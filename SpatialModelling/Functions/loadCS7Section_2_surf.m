function [Dexp1,Arc,Xtest,Output] = loadCS7Section_2_surf(D3,Dexp)

fid = fopen('CS7/Proj_shots_BLEND_high.csv', 'rt');
D0_ = textscan(fid,'%f %f %f %s %s','headerLines', 1,'Delimiter',',');
fclose(fid);
D1.data = cell2mat(D0_(:,1:3));
D1.textdata =D0_{:,4};

Arc = importdata('Data/EmDiscCS7_Arc.mat');
%OBJ3=importdata('Data/CS7.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ3=importdata('Data/CS7_section_cs2.mat');

%Load in the medium quality plots
%D0_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/E15B/Proj_shots_BLEND_high_CP.csv');

%D1_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/EmDisc_BLEND_points_high_final.csv');
%D2_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/EmDisc_BLEND_faces_high_final.csv');

%D3_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/Tb_BLEND_faces_high_final.csv');
%D4_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/Am_BLEND_faces_high_final.csv');
%D5_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/sys_BLEND_faces_high_final.csv'); %Trophoblast
%D6_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/ExMes_BLEND_faces_high_final.csv');
%D7_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/ExMes_stalk_BLEND_faces_high_final.csv');



%innds1 = (find(D6_(:,1)<size(D1_,1)));
%innds2 = (find(D6_(:,2)<size(D1_,1)));
%innds3 = (find(D6_(:,3)<size(D1_,1)));
%Xex = [D1_(D6_(innds1,1),1),D1_(D6_(innds1,1),2),D1_(D6_(innds1,1),3); 
%    D1_(D6_(innds2,2),1),D1_(D6_(innds2,2),2),D1_(D6_(innds2,2),3) ; 
%    D1_(D6_(innds3,3),1),D1_(D6_(innds3,3),2),D1_(D6_(innds3,3),3)];
Xex2 = [OBJ3.vertices(OBJ3.objects(24).data.vertices,:)];


%Xem = [D1_(D2_(:,1),1),D1_(D2_(:,1),2),D1_(D2_(:,1),3); 
%    D1_(D2_(:,2),1),D1_(D2_(:,2),2),D1_(D2_(:,2),3) ; 
 %   D1_(D2_(:,3),1),D1_(D2_(:,3),2),D1_(D2_(:,3),3)];


 
Xem2 = OBJ3.vertices(OBJ3.objects(8).data.vertices,:);


%Xem = Xem(find(Xem(:,2)<600),:);

%Xtroph = [D1_(D3_(:,1),1),D1_(D3_(:,1),2),D1_(D3_(:,1),3); D1_(D3_(:,2),1),D1_(D3_(:,2),2),D1_(D3_(:,2),3) ; D1_(D3_(:,3),1),D1_(D3_(:,3),2),D1_(D3_(:,3),3)];

Xtroph2 = OBJ3.vertices(OBJ3.objects(16).data.vertices,:);


%Xtroph = Xtroph(find(Xtroph(:,1)>-600),:);

%Xam = [D1_(D4_(:,1),1),D1_(D4_(:,1),2),D1_(D4_(:,1),3); D1_(D4_(:,2),1),D1_(D4_(:,2),2),D1_(D4_(:,2),3) ; D1_(D4_(:,3),1),D1_(D4_(:,3),2),D1_(D4_(:,3),3)];

Xam2 = OBJ3.vertices(OBJ3.objects(4).data.vertices,:);


%Xsys = [D1_(D5_(:,1),1),D1_(D5_(:,1),2),D1_(D5_(:,1),3); D1_(D5_(:,2),1),D1_(D5_(:,2),2),D1_(D5_(:,2),3) ; D1_(D5_(:,3),1),D1_(D5_(:,3),2),D1_(D5_(:,3),3)];
Xsys2 = [OBJ3.vertices(OBJ3.objects(12).data.vertices,:)];
   % OBJ3.vertices(OBJ3.objects(12).data.vertices,:)];
%Xexst = [D1_(D7_(:,1),1),D1_(D7_(:,1),2),D1_(D7_(:,1),3); 
%    D1_(D7_(:,2),1),D1_(D7_(:,2),2),D1_(D7_(:,2),3) ; 
%    D1_(D7_(:,3),1),D1_(D7_(:,3),2),D1_(D7_(:,3),3)];

%Xexst2 = OBJ3.vertices(OBJ3.objects(30).data.vertices,:);


%Xve2 = OBJ3.vertices(OBJ3.objects(12).data.vertices,:);


%Arc1 = [ones(1,1000)',-18*ones(1,1000)',linspace(-115,217,1000)'];

%[R,xcyc] = fit_circle_through_3_points([-135 -28; 1 10; 124 -33]);

%thet1 = linspace(0.30*pi,0.66*pi,1000);
%thet2 = linspace(0.30*pi,0.66*pi,1000);
%thet3 = linspace(0.40*pi,0.61*pi,1000);
%Arc2 = [R*cos(thet1)' + xcyc(1), R*sin(thet1)' + xcyc(2), 160*ones(1,1000)'];
%Arc3 = [R*cos(thet2)' + xcyc(1), R*sin(thet2)' + xcyc(2), 40*ones(1,1000)'];
%Arc4 = [R*cos(thet3)' + xcyc(1), R*sin(thet3)' + xcyc(2), -60*ones(1,1000)'];

% plot3(Xsys(:,1),Xsys(:,2),Xsys(:,3),'r.') %Sys
% hold on
% plot3(Xem(:,1),Xem(:,2),Xem(:,3),'b.') %EmDisc
% plot3(Arc1(:,1),Arc1(:,2),Arc1(:,3),'c.-') %EmDisc
% plot3(Arc2(:,1),Arc2(:,2),Arc2(:,3),'c.-') %EmDisc
% plot3(Arc3(:,1),Arc3(:,2),Arc3(:,3),'c.-') %EmDisc
% plot3(Arc4(:,1),Arc4(:,2),Arc4(:,3),'c.-') %EmDisc

%Arc = [Arc1];

%Xex2 = unique(Xex, 'rows');
%Xexst2 = unique(Xexst, 'rows');

%Xam2 = unique(Xam, 'rows');
%Xsys2 = unique(Xsys, 'rows');
%Xem2 = unique(Xem, 'rows');
%Xtroph2 = unique(Xtroph, 'rows');

D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ;  Xtroph2,5*ones(size(Xtroph2,1),1) ];


locations = D1.textdata(1:end,1);
for i = 1:length(locations)
    try
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , locations{i}));   
    %catch
    %idx(i,1) = find(strcmp(D3.textdata(:,5) , ['E25A_',strrep(locations{i},'-','_')]));    
    end    
end
CellIDs = D3.textdata(idx+1,1);
D3subs = D3.data(idx,:);
CellIDsu = D3.textdata(idx+1,6)
CellIDs3 = D3.textdata(idx+1,9);
CellIDs2 = D3.textdata(idx+1,11);
CellIDs4 = D3.textdata(idx+1,5);

%Need to change - to . in CellIDs or vi 

genes = Dexp.textdata(2:end,1);
headers = Dexp.textdata(1,2:end);
headers = strrep(headers,'"','');
for i = 1:length(CellIDs)
    try
    indx(i,1) = find(strcmp(headers , CellIDs{i} ));
    catch
        try
    indx(i,1) = find(strcmp(headers , strrep(CellIDs{i}, '-','.' )));   
        catch
        end
    end
end

[CellIDs(find(indx==0)),CellIDs2(find(indx==0)),CellIDs4(find(indx==0))]

CellIDs = CellIDs(find(indx~=0));
CellIDs2 = CellIDs2(find(indx~=0));
CellIDs3 = CellIDs3(find(indx~=0));
CellIDsu = CellIDsu(find(indx~=0));
CellIDs4 = CellIDs4(find(indx~=0));

XYZ = D1.data(find(indx~=0),:);

Dexp1 = Dexp.data(:,indx(find(indx~=0)));
genes = Dexp.textdata(2:end,1);
Xtrain = XYZ;
Xtrain(:,1:3) = Xtrain(:,1:3) / 400;
Xtest = D2(:,[1:3]);
Xtest(:,1:3) =Xtest(:,1:3)/400;

%Labs = CellIDs3;%D3.textdata(inds,10);
LabBin = NaN*zeros(length(CellIDs2),1);


%Training labels
%keyboard
LabBin(find(strcmp(CellIDs2,'Am_CS7')==1)) = 0;
LabBin(find(strcmp(CellIDs2,'EmDisc_CS7')==1)) = 1;
LabBin(find(strcmp(CellIDs2,'SYS_CS7')==1)) = 2;
LabBin(find(strcmp(CellIDs2,'Tb_CS7')==1)) = 5;
LabBin(find(strcmp(CellIDs2,'ExMes_CS7')==1)) = 3;

ind0 = find(D2(:,4)==0);
ind1 = find(D2(:,4)==1);
ind2 = find(D2(:,4)==2);
ind3 = find(D2(:,4)==3);
ind5 = find(D2(:,4)==5);

Output.Xtest = Xtest;
Output.Xtrain = Xtrain;
Output.LabBin = LabBin;
Output.genes = genes;
Output.ind0 = ind0;
Output.ind1 = ind1;
Output.ind2 = ind2;
Output.ind3 = ind3;
Output.ind5 = ind5;
Output.CellIDs = CellIDs4;

return

% 
% %keyboard
%keyboard
alpha = 0.05;
 DMAP1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/CCAFigForESCs/DimRed/ESC.csv')
 estarg = DMAP1.textdata(1,2:end);
 estarg = strrep(estarg,'"','');
 
 DMAP1.textdata(2:end,1);
 targ = DMAP1.textdata(2:end,1);
 for ii = 1:length(estarg)
 Targ = targ{find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii)))};
 Targ = strrep(Targ,'EmDisc_CS7_','');
 Targ = strrep(Targ,'ExMes_CS7_','');
 Targ = strrep(Targ,'VE_CS7_','');
 Targ = strrep(Targ,'Tb_CS7_','');
 Targ = strrep(Targ,'Am_CS7_','');
 try
 Ind0(ii) = find(strcmp(CellIDs4,Targ)==1)
 P0(ii) = DMAP1.data(find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii))),ii);
 catch
 end
 end
 
 try
 P01 = P0(Ind0~=0);
 UMapped = unique(CellIDs4(Ind0(Ind0~=0)));
 
 indss = [];
 for i = 1:length(UMapped)
 i0 = find(strcmp(CellIDs4,UMapped{i})==1);
 indss = [indss; i0];
 end
 
 

 subplot(4,12,[1 2 3 4 13 14 15 16 25 26 27 28 37 38 39 40]);
 h = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill', 'MarkerEdgeColor',[2,51,191]/255, 'MarkerFaceColor',[2,51,191]/255)
 hold on
 scatter3(Output.Xtrain(setdiff(find(Output.LabBin==1),indss),1),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),2),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),3),20,'fill','MarkerEdgeColor','k', 'MarkerFaceColor','k')

% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),20,'fill')
 %text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
 %scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
 %scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
 %text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
 scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill','MarkerEdgeColor','r', 'MarkerFaceColor','r')
% text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([11.3486,74.1824])
 set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
 set(h,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)

 catch
     subplot(4,12,[1 2 3 4 13 14 15 16 25 26 27 28 37 38 39 40]);
 h = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill', 'MarkerEdgeColor',[2,51,191]/255, 'MarkerFaceColor',[2,51,191]/255)
 hold on
 scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),20,'fill', 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
 %text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
 %scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
 %scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
 %text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
 %scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill')
 %text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([11.3486,74.1824])
 set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
 set(h,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)

 end
 
 
 DMAP1.textdata(2:end,1);
 Ind0 = [];
 targ = DMAP1.textdata(2:end,1);
 for ii = 1:length(estarg)
 Targ = targ(find(DMAP1.data(:,ii)>0.06));
 try
     for kk = 1:length(Targ)
         try
             Targ1 = strrep(Targ{kk},'EmDisc_CS7_','');
 Targ1 = strrep(Targ1,'ExMes_CS7_','');
 Targ1 = strrep(Targ1,'SYS_CS7_','');
 Targ1 = strrep(Targ1,'Tb_CS7_','');
 Targ1 = strrep(Targ1,'Am_CS7_','');
     Ind0 = [Ind0;find(strcmp(CellIDs4,Targ1)==1)];
 % = DMAP1.data(find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii))),ii);
         catch
         end
     end
 end
 end
 

 
 try
 
 P01 = P0(Ind0~=0);
 UMapped = unique(CellIDs4(Ind0(Ind0~=0)));
 
 indss = [];
 for i = 1:length(UMapped)
 i0 = find(strcmp(CellIDs4,UMapped{i})==1);
 indss = [indss; i0];
 end
 
 
 
 subplot(4,12,[5 6 7 8 17 18 19 20 29 30 31 32 41 42 43 44]);
 h2 = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill', 'MarkerEdgeColor',[2,51,191]/255, 'MarkerFaceColor',[2,51,191]/255)
 hold on
  scatter3(Output.Xtrain(setdiff(find(Output.LabBin==1),indss),1),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),2),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),3),20,'fill','MarkerEdgeColor','k', 'MarkerFaceColor','k')

% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),20,'fill', 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
 %text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
 %scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
 %scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
 %text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
 scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill', 'MarkerEdgeColor','r', 'MarkerFaceColor','r')
 %text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([11.3486,74.1824])
 set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
 
set(h2,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)

 catch

 indss = [];
 subplot(4,12,[5 6 7 8 17 18 19 20 29 30 31 32 41 42 43 44]);
 h2 = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill', 'MarkerEdgeColor',[2,51,191]/255, 'MarkerFaceColor',[2,51,191]/255)
 hold on
  scatter3(Output.Xtrain(setdiff(find(Output.LabBin==1),indss),1),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),2),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),3),20,'fill','MarkerEdgeColor','k', 'MarkerFaceColor','k')

% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),20,'fill', 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
 %text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
 %scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
 %scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
 %text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
% scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill', 'MarkerEdgeColor','r', 'MarkerFaceColor','r')
 %text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([11.3486,74.1824])
 set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
 
set(h2,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
     
     
 end

% 
% h = figure(1)
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
% savefig(h,['~/Desktop/FIG/CS7_EmDisc'])
% clf
% 
% h = figure(2)
% scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
% savefig(h,['~/Desktop/FIG/CS7_Am'])
% clf
% 
% h = figure(3)
% scatter3(Output.Xtest(Output.ind2,1),Output.Xtest(Output.ind2,2),Output.Xtest(Output.ind2,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==2),1),Output.Xtrain(find(Output.LabBin==2),2),Output.Xtrain(find(Output.LabBin==2),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==2),1),Output.Xtrain(find(Output.LabBin==2),2),Output.Xtrain(find(Output.LabBin==2),3), CellIDs4(find(Output.LabBin==2)));
% savefig(h,['~/Desktop/FIG/CS7_SYS'])
% clf
% 
% h = figure(4)
% scatter3(Output.Xtest(Output.ind3,1),Output.Xtest(Output.ind3,2),Output.Xtest(Output.ind3,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3), CellIDs4(find(Output.LabBin==3)));
% savefig(h,['~/Desktop/FIG/CS7_ExMes'])
% clf
% 
% h = figure(5)
% scatter3(Output.Xtest(Output.ind5,1),Output.Xtest(Output.ind5,2),Output.Xtest(Output.ind5,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3), CellIDs4(find(Output.LabBin==5)));
% savefig(h,['~/Desktop/FIG/CS7_Tb'])
% clf
% 
% 
% 
% 
% h = figure(6)
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
% scatter3(Output.Xtest(Output.ind3,1),Output.Xtest(Output.ind3,2),Output.Xtest(Output.ind3,3),5,'fill')
% scatter3(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3), CellIDs4(find(Output.LabBin==3)));
% 
% savefig(h,['~/Desktop/FIG/CS7_EmDiscExMes'])
% %clf

%figure(4); 
%scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
%hold on
%scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')




%xtest = [Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3)]


%CellIDs(find(Output.LabBin==1),1)
%CellIDs2(find(Output.LabBin==1),1)
