function [Dexp1,Arc,Xtest,Output] = loadCS6(D3,Dexp)


fid = fopen('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/Proj_shots_BLEND_high_update.csv', 'rt'); 
D0_ = textscan(fid,'%f %f %f %s %s','headerLines', 1,'Delimiter',',');
fclose(fid);
D1.data = cell2mat(D0_(:,1:3));
D1.textdata =D0_{:,4};

%Load in the medium quality plots
%D0_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/E15B/Proj_shots_BLEND_high_CP.csv');

D1_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/EmDisc_BLEND_points_high_final.csv');
D2_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/EmDisc_BLEND_faces_high_final.csv');
D3_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/Tb_BLEND_faces_high_final.csv');
D4_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/Am_BLEND_faces_high_final.csv');
D5_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/sys_BLEND_faces_high_final.csv'); %Trophoblast
D6_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/ExMes_BLEND_faces_high_final.csv');
D7_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/VE_BLEND_faces_high_final.csv');
D8_ = importdata('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/PGC_BLEND_faces_high_final.csv');



Xex = [D1_(D6_(:,1),1),D1_(D6_(:,1),2),D1_(D6_(:,1),3); D1_(D6_(:,2),1),D1_(D6_(:,2),2),D1_(D6_(:,2),3) ; D1_(D6_(:,3),1),D1_(D6_(:,3),2),D1_(D6_(:,3),3)];

Xem = [D1_(D2_(:,1),1),D1_(D2_(:,1),2),D1_(D2_(:,1),3); D1_(D2_(:,2),1),D1_(D2_(:,2),2),D1_(D2_(:,2),3) ; D1_(D2_(:,3),1),D1_(D2_(:,3),2),D1_(D2_(:,3),3)];

Xtroph = [D1_(D3_(:,1),1),D1_(D3_(:,1),2),D1_(D3_(:,1),3); D1_(D3_(:,2),1),D1_(D3_(:,2),2),D1_(D3_(:,2),3) ; D1_(D3_(:,3),1),D1_(D3_(:,3),2),D1_(D3_(:,3),3)];



Xam = [D1_(D4_(:,1),1),D1_(D4_(:,1),2),D1_(D4_(:,1),3); D1_(D4_(:,2),1),D1_(D4_(:,2),2),D1_(D4_(:,2),3) ; D1_(D4_(:,3),1),D1_(D4_(:,3),2),D1_(D4_(:,3),3)];

%Somethinng weird with SYS & VE?
Xsys = [D1_(D5_(:,1),1),D1_(D5_(:,1),2),D1_(D5_(:,1),3); D1_(D5_(:,2),1),D1_(D5_(:,2),2),D1_(D5_(:,2),3) ; D1_(D5_(:,3),1),D1_(D5_(:,3),2),D1_(D5_(:,3),3)];
Xsys = Xsys(find(Xsys(:,2)>-400),:);
Xsys = Xsys(find( abs(Xsys(:,1)-218.5)>1 & abs(Xsys(:,2)-128.4)>1),:);

Xve = [D1_(D7_(:,1),1),D1_(D7_(:,1),2),D1_(D7_(:,1),3); D1_(D7_(:,2),1),D1_(D7_(:,2),2),D1_(D7_(:,2),3) ; D1_(D7_(:,3),1),D1_(D7_(:,3),2),D1_(D7_(:,3),3)];

Xve = Xve(find(Xve(:,2)>-400),:);
Xve = Xve(find(Xve(:,1)<350),:);

Xpgc = [D1_(D8_(:,1),1),D1_(D8_(:,1),2),D1_(D8_(:,1),3); D1_(D8_(:,2),1),D1_(D8_(:,2),2),D1_(D8_(:,2),3) ; D1_(D8_(:,3),1),D1_(D8_(:,3),2),D1_(D8_(:,3),3)];

Arc1 = [ones(1,1000)',-18*ones(1,1000)',linspace(-115,217,1000)'];

[R,xcyc] = fit_circle_through_3_points([-135 -28; 1 10; 124 -33])

%thet1 = linspace(0.30*pi,0.66*pi,1000);
%thet2 = linspace(0.30*pi,0.66*pi,1000);
%thet3 = linspace(0.40*pi,0.61*pi,1000);
%Arc2 = [R*cos(thet1)' + xcyc(1), R*sin(thet1)' + xcyc(2), 160*ones(1,1000)'];
%Arc3 = [R*cos(thet2)' + xcyc(1), R*sin(thet2)' + xcyc(2), 40*ones(1,1000)'];
%Arc4 = [R*cos(thet3)' + xcyc(1), R*sin(thet3)' + xcyc(2), -60*ones(1,1000)'];

Arc = [Arc1];

Xex2 = unique(Xex, 'rows');
Xam2 = unique(Xam, 'rows');
Xve2 = unique(Xve, 'rows');
Xsys2 = unique(Xsys, 'rows');
Xem2 = unique(Xem, 'rows');
Xtroph2 = unique(Xtroph, 'rows');
Xpgc2 = unique(Xpgc, 'rows');


D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ; Xve2,4*ones(size(Xve2,1),1) ; Xtroph2,5*ones(size(Xtroph2,1),1); Xpgc2,6*ones(size(Xpgc2,1),1)  ];

locations = D1.textdata(1:end,1)
for i = 1:length(locations)
    try
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , ['E15C_',strrep(locations{i},'-','_')]));   
    end    
end
CellIDs = D3.textdata(idx+1,1);
D3subs = D3.data(idx,:);
CellIDsu = D3.textdata(idx+1,6)
CellIDs2 = D3.textdata(idx+1,11);
CellIDs4 = D3.textdata(idx+1,5);

%CellIDs = D3.textdata(idx+1,1);
%CellIDsu = D3.textdata(idx+1,6)
%D3subs = D3.data(idx,:);
%CellIDs2 = D3.textdata(idx+1,11);
%CellIDs4 = D3.textdata(idx+1,5);



%Load in expression data
Dexp = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Final_AllGoodShots_wCS6r/NormData.csv');
genes = Dexp.textdata(2:end,1);
headers = Dexp.textdata(1,2:end);
headers = strrep(headers,'"','');
headers = strrep(headers,'X3536STDY','3536STDY');

for i = 1:length(CellIDs)
    try
    indx(i,1) = find(strcmp(headers , CellIDs{i} ));
    catch
    end
end



%Are there any that don't exist?
CellIDs = CellIDs(find(indx~=0));
CellIDs2 = CellIDs2(find(indx~=0))
CellIDs4 = CellIDs4(find(indx~=0))
CellIDsu = CellIDsu(find(indx~=0))

XYZ = D1.data(find(indx~=0),:);


Dexp1 = Dexp.data(:,indx(find(indx~=0)));
genes = Dexp.textdata(2:end,1);
Xtrain = XYZ;
Xtrain(:,1:3) = Xtrain(:,1:3) / 400;
Xtest = D2(:,[1:3]);
Xtest(:,1:3) =Xtest(:,1:3)/400;
%Xtest2 = D2(:,[1:3]);
%Xtest2(:,1:3) =Xtest2(:,1:3)/400;


%Now get ...
%for ii = 1:length(CellIDs)
%      inds(ii,1) = find(strcmp(CellIDs{ii},D3.textdata(1:end,1))==1);    
%end

%Labs = D3.textdata(inds,10);
LabBin = NaN*zeros(length(CellIDs2),1);

%Training labels
LabBin(find(strcmp(CellIDs2,'Am_CS6')==1)) = 0;
LabBin(find(strcmp(CellIDs2,'EmDisc_CS6')==1)) = 1;
LabBin(find(strcmp(CellIDs2,'VE_CS6')==1)) = 4;
LabBin(find(strcmp(CellIDs2,'SYS_CS6')==1)) = 2;
LabBin(find(strcmp(CellIDs2,'Tb_CS6')==1)) = 5;
LabBin(find(strcmp(CellIDs2,'ExMes_CS6')==1)) = 3;
LabBin(find(strcmp(CellIDs2,'PGC_CS6')==1)) = 6;

%Loop and get cell iDs for tissue specific models

%D2 = [Xprot2,3*ones(size(Xprot2,1),1) ; Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xex2,2*ones(size(Xex2,1),1) ; Xve2,3*ones(size(Xve2,1),1) ; Xtroph2,4*ones(size(Xtroph2,1),1) ];


%D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ; Xve2,4*ones(size(Xve2,1),1) ; Xtroph2,5*ones(size(Xtroph2,1),1) ];
%keyboard
ind0 = find(D2(:,4)==0);
ind1 = find(D2(:,4)==1);
ind2 = find(D2(:,4)==2);
ind3 = find(D2(:,4)==3);
ind4 = find(D2(:,4)==4);
ind5 = find(D2(:,4)==5);
ind6 = find(D2(:,4)==6);

Output.Xtest = Xtest;
Output.Xtrain = Xtrain;
Output.LabBin = LabBin;
Output.genes = genes;
Output.ind0 = ind0;
Output.ind1 = ind1;
Output.ind2 = ind2;
Output.ind3 = ind3;
Output.ind4 = ind4;
Output.ind5 = ind5;
Output.ind6 = ind6;
return
%keyboard
DMAP1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/CCAFigForESCs/DimRed/ESC_pi.csv')
estarg = DMAP1.textdata(1,2:end);
estarg = strrep(estarg,'"','');

DMAP1.textdata(2:end,1);
targ = DMAP1.textdata(2:end,1);
for ii = 1:length(estarg)
Targ = targ{find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii)))};
Targ = strrep(Targ,'EmDisc_CS6_','');
Targ = strrep(Targ,'ExMes_CS6_','');
Targ = strrep(Targ,'VE_CS6_','');
Targ = strrep(Targ,'Tb_CS6_','');
Targ = strrep(Targ,'Am_CS6_','');
Targ = strrep(Targ,'PGC_CS6_','');

try
Ind0(ii) = find(strcmp(CellIDs4,Targ)==1)
P0(ii) = DMAP1.data(find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii))),ii);
catch
end
end

P01 = P0(Ind0~=0);
UMapped = unique(CellIDs4(Ind0(Ind0~=0)));

indss = [];
for i = 1:length(UMapped)
i0 = find(strcmp(CellIDs4,UMapped{i})==1);
indss = [indss; i0];
end

subplot(4,12,[1 2 3 4 13 14 15 16 25 26 27 28 37 38 39 40]);
scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
hold on
scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
%scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
%scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill')
text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([12.8199,21.3768])
set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))

DMAP1.textdata(2:end,1);
Ind0 = [];
targ = DMAP1.textdata(2:end,1);
for ii = 1:length(estarg)
Targ = targ(find(DMAP1.data(:,ii)>0.06));
try
    for kk = 1:length(Targ)
        try
            Targ1 = strrep(Targ{kk},'EmDisc_CS6_','');
Targ1 = strrep(Targ1,'ExMes_CS6_','');
Targ1 = strrep(Targ1,'VE_CS6_','');
Targ1 = strrep(Targ1,'Tb_CS6_','');
Targ1 = strrep(Targ1,'Am_CS6_','');
Targ1 = strrep(Targ1,'PGC_CS6_','');
    Ind0 = [Ind0;find(strcmp(CellIDs4,Targ1)==1)];
% = DMAP1.data(find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii))),ii);
        catch
        end
    end
end
end
%P01 = P0(Ind0~=0);
UMapped = unique(CellIDs4(Ind0(Ind0~=0)));

indss = [];
for i = 1:length(UMapped)
i0 = find(strcmp(CellIDs4,UMapped{i})==1);
indss = [indss; i0];
end



subplot(4,12,[5 6 7 8 17 18 19 20 29 30 31 32 41 42 43 44]);
scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
hold on
scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
%scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
%scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill')
text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([12.8199,21.3768])
set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))



Output.CellIDs = CellIDs4;

% 
% h = figure(1)
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
% scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
% scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));

%keyboard

%savefig(h,['~/Desktop/FIG/CS6_Am'])
%clf

% 
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
% hold on
% scatter3(Output.Xtest(Output.ind6,1),Output.Xtest(Output.ind6,2),Output.Xtest(Output.ind6,3),5,'fill')
% 
% 
% 
% h = figure(1)
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
% savefig(h,['~/Desktop/FIG/CS6_EmDisc'])
% clf
% 
% 
% 
% 
% 
% h = figure(2)
% scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
% savefig(h,['~/Desktop/FIG/CS6_Am'])
% clf
% 
% h = figure(3)
% scatter3(Output.Xtest(Output.ind2,1),Output.Xtest(Output.ind2,2),Output.Xtest(Output.ind2,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==2),1),Output.Xtrain(find(Output.LabBin==2),2),Output.Xtrain(find(Output.LabBin==2),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==2),1),Output.Xtrain(find(Output.LabBin==2),2),Output.Xtrain(find(Output.LabBin==2),3), CellIDs4(find(Output.LabBin==2)));
% savefig(h,['~/Desktop/FIG/CS6_SYS'])
% clf
% 
% h = figure(4)
% scatter3(Output.Xtest(Output.ind3,1),Output.Xtest(Output.ind3,2),Output.Xtest(Output.ind3,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3), CellIDs4(find(Output.LabBin==3)));
% savefig(h,['~/Desktop/FIG/CS6_ExMes'])
% clf
% 
% h = figure(4)
% scatter3(Output.Xtest(Output.ind4,1),Output.Xtest(Output.ind4,2),Output.Xtest(Output.ind4,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==4),1),Output.Xtrain(find(Output.LabBin==4),2),Output.Xtrain(find(Output.LabBin==4),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==4),1),Output.Xtrain(find(Output.LabBin==4),2),Output.Xtrain(find(Output.LabBin==4),3), CellIDs4(find(Output.LabBin==4)));
% savefig(h,['~/Desktop/FIG/CS6_VE'])
% clf
% 
% h = figure(5)
% scatter3(Output.Xtest(Output.ind5,1),Output.Xtest(Output.ind5,2),Output.Xtest(Output.ind5,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3), CellIDs4(find(Output.LabBin==5)));
% savefig(h,['~/Desktop/FIG/CS6_Tb'])
% clf
% 
% 
% h = figure(6)
% scatter3(Output.Xtest(Output.ind6,1),Output.Xtest(Output.ind6,2),Output.Xtest(Output.ind6,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==6),1),Output.Xtrain(find(Output.LabBin==6),2),Output.Xtrain(find(Output.LabBin==6),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==6),1),Output.Xtrain(find(Output.LabBin==6),2),Output.Xtrain(find(Output.LabBin==6),3), CellIDs4(find(Output.LabBin==6)));
% savefig(h,['~/Desktop/FIG/CS6_PGC'])
% clf


