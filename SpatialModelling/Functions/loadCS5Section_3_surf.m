function [Dexp1,Arc,Xtest,Output] = loadCS5Section_3_surf(D3,Dexp)

fid = fopen('CS5/Proj_shots_BLEND_high.csv', 'rt'); 
D0_ = textscan(fid,'%f %f %f %s %s','headerLines', 1,'Delimiter',',');
fclose(fid);
D1.data = cell2mat(D0_(:,1:3));
D1.textdata =D0_{:,4};

Arc = importdata('Data/EmDiscCS5_Arc.mat');
OBJ1=importdata('Data/CS5_section_cs4.mat');

Xex2 = OBJ1.vertices(OBJ1.objects(12).data.vertices,:);
Xem2 = OBJ1.vertices(OBJ1.objects(8).data.vertices,:);
Xtroph2 = OBJ1.vertices(OBJ1.objects(20).data.vertices,:);
Xam2 = OBJ1.vertices(OBJ1.objects(4).data.vertices,:);
Xsys2 = OBJ1.vertices(OBJ1.objects(16).data.vertices,:);
Xve2 = OBJ1.vertices(OBJ1.objects(24).data.vertices,:);
D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ; Xve2,4*ones(size(Xve2,1),1) ; Xtroph2,5*ones(size(Xtroph2,1),1) ];


% D2h = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xex2,2*ones(size(Xex2,1),1) ; Xve2,3*ones(size(Xve2,1),1) ; Xtroph2,4*ones(size(Xtroph2,1),1) ];
locations = D1.textdata(1:end,1)
for i = 1:length(locations)
    try
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , ['P1_E15B_',strrep(locations{i},'-','_')]));   
    catch
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , ['P2_E15B_',strrep(locations{i},'-','_')]));    
    end    
end
CellIDs = D3.textdata(idx+1,1);
CellIDsu = D3.textdata(idx+1,6)
D3subs = D3.data(idx,:);
CellIDs2 = D3.textdata(idx+1,11);
CellIDs4 = D3.textdata(idx+1,5);


%Load in expression data
genes = Dexp.textdata(2:end,1);
headers = Dexp.textdata(1,2:end);
headers = strrep(headers,'"','');
for i = 1:length(CellIDs)
    try
    indx(i,1) = find(strcmp(headers , CellIDs{i} ));
    catch
    end
end

CellIDs = CellIDs(find(indx~=0));
CellIDs2 = CellIDs2(find(indx~=0));
CellIDs4 = CellIDs4(find(indx~=0));
CellIDsu = CellIDsu(find(indx~=0))

XYZ = D1.data(find(indx~=0),:);

Dexp1 = Dexp.data(:,indx(find(indx~=0)));
genes = Dexp.textdata(2:end,1);
Xtrain = XYZ;
Xtrain(:,1:3) = Xtrain(:,1:3) / 400;
Xtest = D2(:,[1:3]);
Xtest(:,1:3) =Xtest(:,1:3)/400;

LabBin = NaN*zeros(length(CellIDs2),1);

%Training labels
LabBin(find(strcmp(CellIDs2,'Am_CS5')==1)) = 0;
LabBin(find(strcmp(CellIDs2,'EmDisc_CS5')==1)) = 1;
LabBin(find(strcmp(CellIDs2,'VE_CS5')==1)) = 4;
LabBin(find(strcmp(CellIDs2,'SYS_CS5')==1)) = 2;
LabBin(find(strcmp(CellIDs2,'Tb_CS5')==1)) = 5;
LabBin(find(strcmp(CellIDs2,'ExMes_CS5')==1)) = 3;

ind0 = find(D2(:,4)==0);
ind1 = find(D2(:,4)==1);
ind2 = find(D2(:,4)==2);
ind3 = find(D2(:,4)==3);
ind4 = find(D2(:,4)==4);
ind5 = find(D2(:,4)==5);

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

Output.CellIDs = CellIDs4;

return

DMAP1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/CCAFigForESCs/DimRed/ESC.csv')
estarg = DMAP1.textdata(1,2:end);
estarg = strrep(estarg,'"','');

DMAP1.textdata(2:end,1);
targ = DMAP1.textdata(2:end,1);
for ii = 1:length(estarg)
Targ = targ{find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii)))};
Targ = strrep(Targ,'EmDisc_CS5_','');
Targ = strrep(Targ,'ExMes_CS5_','');
Targ = strrep(Targ,'VE_CS5_','');
Targ = strrep(Targ,'Tb_CS5_','');
Targ = strrep(Targ,'Am_CS5_','');
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
h = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill','MarkerEdgeColor',[12,156,245]/255, 'MarkerFaceColor',[12,156,245]/255)
hold on
%scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')

%text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
%scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
scatter3(Output.Xtrain(setdiff(find(Output.LabBin==1),indss),1),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),2),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),3),20,'fill','MarkerEdgeColor','k', 'MarkerFaceColor','k')
%text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill','MarkerEdgeColor','r', 'MarkerFaceColor','r')
%text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([ -79.1408,77.6211])
%view([-78.2671,26.9624])
set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
alpha = 0.05;
set(h,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
%set(h2,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
catch
    subplot(4,12,[1 2 3 4 13 14 15 16 25 26 27 28 37 38 39 40]);
h = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill','MarkerEdgeColor',[12,156,245]/255, 'MarkerFaceColor',[12,156,245]/255)
hold on
%scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
indss = [];
%text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
%scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
scatter3(Output.Xtrain(setdiff(find(Output.LabBin==1),indss),1),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),2),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),3),20,'fill','MarkerEdgeColor','k', 'MarkerFaceColor','k')
%text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
%scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill','MarkerEdgeColor','r', 'MarkerFaceColor','r')
%text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
view([ -79.1408,77.6211])
%view([-78.2671,26.9624])
set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
alpha = 0.05;
set(h,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
%set(h2,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
end


DMAP1.textdata(2:end,1);
Ind0 = [];
targ = DMAP1.textdata(2:end,1);
for ii = 1:length(estarg)
Targ = targ(find(DMAP1.data(:,ii)>0.06));
try
    for kk = 1:length(Targ)
        try
            Targ1 = strrep(Targ{kk},'EmDisc_CS5_','');
Targ1 = strrep(Targ1,'ExMes_CS5_','');
Targ1 = strrep(Targ1,'VE_CS5_','');
Targ1 = strrep(Targ1,'Tb_CS5_','');
Targ1 = strrep(Targ1,'Am_CS5_','');
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

%subplot(4,12,[5 6 7 8 17 18 19 20 29 30 31 32 41 42 43 44]);
%h2 = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill','MarkerEdgeColor',[12,156,245]/255, 'MarkerFaceColor',[12,156,245]/255)
%hold on
%scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
%scatter3(Output.Xtrain(setdiff(Output.Xtrain(indss,1),find(Output.LabBin==0),1)),Output.Xtrain(setdiff(Output.Xtrain(indss,1),find(Output.LabBin==0)),2),Output.Xtrain(Output.Xtrain(indss,1),setdiff(find(Output.LabBin==0)),3),20,'fill')
%scatter3(Output.Xtrain(setdiff(find(Output.LabBin==1),indss),1),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),2),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),3),20,'fill','MarkerEdgeColor','k', 'MarkerFaceColor','k')

%text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
%scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
%scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
%scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill','MarkerEdgeColor','r', 'MarkerFaceColor','r')
%text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
%view([ -79.1408,77.6211])
%set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
%set(h2,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)

% keyboard
% 
% 
% 
% 
% DMAP1.textdata(2:end,1);
% 
% targ = DMAP1.textdata(2:end,1);
% for ii = 1:length(estarg)
% Targ = targ(find(DMAP1.data(:,ii)>0.06));
% Ind0 = [];
% try
%     for kk = 1:length(Targ)
%         try
%             Targ1 = strrep(Targ{kk},'EmDisc_CS5_','');
% Targ1 = strrep(Targ1,'ExMes_CS5_','');
% Targ1 = strrep(Targ1,'VE_CS5_','');
% Targ1 = strrep(Targ1,'Tb_CS5_','');
% Targ1 = strrep(Targ1,'Am_CS5_','');
%     Ind0 = [Ind0;find(strcmp(CellIDs4,Targ1)==1)];
% % = DMAP1.data(find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii))),ii);
%         catch
%         end
%     end
% end
% 
% UMapped = unique(CellIDs4(Ind0(Ind0~=0)));
% indss = [];
% for i = 1:length(UMapped)
% i0 = find(strcmp(CellIDs4,UMapped{i})==1);
% indss = [indss; i0];
% end
% 
% subplot(5,5,ii);
% %subplot(4,12,[5 6 7 8 17 18 19 20 29 30 31 32 41 42 43 44]);
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
% %text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
% %scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
% %scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
% %text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
% scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill')
% text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
% view([-78.2671,26.9624])
% set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% 
% 
% end
%P01 = P0(Ind0~=0);


DMAP1.textdata(2:end,1);
P = zeros(length(CellIDs4),length(estarg));

targ = DMAP1.textdata(2:end,1);
Ind0 = zeros(length(targ),1);
for ii = 1:length(estarg)
Targ = targ(find(DMAP1.data(:,ii)>-1));
    for kk = 1:length(Targ)
        try
Targ1 = strrep(Targ{kk},'EmDisc_CS5_','');
Targ1 = strrep(Targ1,'ExMes_CS5_','');
Targ1 = strrep(Targ1,'VE_CS5_','');
Targ1 = strrep(Targ1,'Tb_CS5_','');
Targ1 = strrep(Targ1,'Am_CS5_','');
Targ1 = strrep(Targ1,'SYS_CS5_','');
    Ind0 = [Ind0;find(strcmp(CellIDs4,Targ1)==1)];
    
    %find(strcmp(Targ1,CellIDs4)==10
    
    %find(strcmp(CellIDs4,targ)==1)
    
    P(find(strcmp(CellIDs4,Targ1)==1),ii) = DMAP1.data(kk,ii);
    
    
    %P0 = DMAP1.data(find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii))),ii);
    
% = DMAP1.data(find(DMAP1.data(:,ii)==max(DMAP1.data(:,ii))),ii);
        catch
        end
    end
end

Output.P = P;

UMapped = targ(find(DMAP1.data(:,ii)>-1));
indss = [];
for i = 1:length(UMapped)
Targ1 = strrep(UMapped{i},'EmDisc_CS5_','');
Targ1 = strrep(Targ1,'ExMes_CS5_','');
Targ1 = strrep(Targ1,'VE_CS5_','');
Targ1 = strrep(Targ1,'Tb_CS5_','');
Targ1 = strrep(Targ1,'Am_CS5_','');
Targ1 = strrep(Targ1,'SYS_CS5_','');    
i0 = find(strcmp(CellIDs4,Targ1)==1);
indss = [indss; i0];
end

XYZP = [Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3)];
Output.XYZP = XYZP;
%UMapped = unique(CellIDs4(Ind0(Ind0~=0)));
%indss = [];
%for i = 1:length(UMapped)
%i0 = find(strcmp(CellIDs4,UMapped{i})==1);
%indss = [indss; i0];
%end

%subplot(5,5,ii);
%subplot(4,12,[5 6 7 8 17 18 19 20 29 30 31 32 41 42 43 44]);
%scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
%hold on
%scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
%scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
%scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
%scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill')
% text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
% view([-78.2671,26.9624])
% set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))


%end

%[1,2,3,4, 5,6,7,8, 9,10,11,12]
%[13,14,15,16,  17,18,19,20,  21,22,23,24]
%[25,26,27,28,  29,30,31,32,  33,34,35,36]
%[37,38,39,40,  41,42,43,44,  45,46,47,48]

%5 6 7 8 17 18 19 20 29 30 31 32 41 42 43 44

%1,2,3,4,13,14,15,16,25,26,27,28,37,38,39,40
%[1,2,3,4, 5,6,7,8, 9,10,11,12]
%[1,2,3,4, 5,6,7,8, 9,10,11,12]

%[1,2,3,4, 5,6,7,8, 9,10,11,12]

%[10,11,12, 13,14,15, 16,17,18]
%[19,20,21, 22,23,24, 25,26,27]

%DMAP2 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/CCAFigForESCs/DimRed/ESC2_pi.csv')

% 
% i0 = find(strcmp(CellIDs4,'P1_E15B_207_H')==1)
% i1 = find(strcmp(CellIDs4,'P1_E15B_207_G')==1)
% i2 = find(strcmp(CellIDs4,'P1_E15B_207_L')==1)
% i3 = find(strcmp(CellIDs4,'P2_E15B_204_F')==1)
% i4 = find(strcmp(CellIDs4,'P1_E15B_207_N')==1) %
% i5 = find(strcmp(CellIDs4,'P1_E15B_207_G')==1)
% i6 = find(strcmp(CellIDs4,'P1_E15B_207_K')==1)
% i7 = find(strcmp(CellIDs4,'P1_E15B_207_C')==1)
% indss = [i1;i2;i3;i4;i5;i6;i7];
% subplot(1,2,2)
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
% scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
% scatter3(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
% scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill')
% 
% 
% % 
% % savefig(h,['~/Desktop/FIG/CS5_Am'])
% % clf
% % 
% % h = figure(3)
% % scatter3(Output.Xtest(Output.ind2,1),Output.Xtest(Output.ind2,2),Output.Xtest(Output.ind2,3),5,'fill')
% % hold on
% % scatter3(Output.Xtrain(find(Output.LabBin==2),1),Output.Xtrain(find(Output.LabBin==2),2),Output.Xtrain(find(Output.LabBin==2),3),95,'fill')
% % text(Output.Xtrain(find(Output.LabBin==2),1),Output.Xtrain(find(Output.LabBin==2),2),Output.Xtrain(find(Output.LabBin==2),3), CellIDs4(find(Output.LabBin==2)));
% % savefig(h,['~/Desktop/FIG/CS5_SYS'])
% % clf
% % 
% % h = figure(4)
% % scatter3(Output.Xtest(Output.ind3,1),Output.Xtest(Output.ind3,2),Output.Xtest(Output.ind3,3),5,'fill')
% % hold on
% % scatter3(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3),95,'fill')
% % text(Output.Xtrain(find(Output.LabBin==3),1),Output.Xtrain(find(Output.LabBin==3),2),Output.Xtrain(find(Output.LabBin==3),3), CellIDs4(find(Output.LabBin==3)));
% % savefig(h,['~/Desktop/FIG/CS5_ExMes'])
% % clf
% % 
% % h = figure(4)
% % scatter3(Output.Xtest(Output.ind4,1),Output.Xtest(Output.ind4,2),Output.Xtest(Output.ind4,3),5,'fill')
% % hold on
% % scatter3(Output.Xtrain(find(Output.LabBin==4),1),Output.Xtrain(find(Output.LabBin==4),2),Output.Xtrain(find(Output.LabBin==4),3),95,'fill')
% % text(Output.Xtrain(find(Output.LabBin==4),1),Output.Xtrain(find(Output.LabBin==4),2),Output.Xtrain(find(Output.LabBin==4),3), CellIDs4(find(Output.LabBin==4)));
% savefig(h,['~/Desktop/FIG/CS5_VE'])
% clf
% 
% h = figure(5)
% scatter3(Output.Xtest(Output.ind5,1),Output.Xtest(Output.ind5,2),Output.Xtest(Output.ind5,3),5,'fill')
% hold on
% scatter3(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3),95,'fill')
% text(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3), CellIDs4(find(Output.LabBin==5)));
% savefig(h,['~/Desktop/FIG/CS5_Tb'])
% clf
% 
