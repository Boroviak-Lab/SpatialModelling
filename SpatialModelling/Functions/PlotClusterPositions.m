Cl = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset_newest_Multiple/Cl9.csv')
Location = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset_newest_Multiple/LOC.csv')

Cl = Cl.data(1:end,1);
Location = Location.textdata(2:end,2);


Location = strrep( Location , 'P1_' , '' );
Location = strrep( Location , 'P2_' , '' );

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

gene = unique({'SOX2','NANOG','POU5F1','MIXL1','T','EOMES','GATA6','LEF1','PDGFRA','ID2','NODAL','CER1','LEFTY2','LHX1','OTX2','HHEX','NANOS3','TFAP2A','TFAP2C','VTCN1','TTR','GC','GATA2','GATA3','HAND2','HGF','ID1','ID3','WNT5A','WT5B','WNT6','TCF4','DKK1','ID1','ID2','NOG','BAMBI','BMP2','BMP6','SOX17','APOA1','GC','ARL13B','IHH','PCH2','SMO','GPC4','GLI1','SOX17','KLF4','DAZL','MAEL','PRAME','NODAL','TDGF1','BMP4','ID2','SOX2','POU5F1','PRDM14','T','LEF1','SFRP1','SFRP2','SFRP5','NANOG','WNT3','WNT8A','FBXO2','TDGF1','MIXL1','GATA6','SOX17','CDH2','CER1','LEFTY2','LHX1','OTX2','NODAL','SDC4','PDGFRA','PDGFA','KDR','VEGFA','DAZL','DDX4','MAEL','DMRT1','OTX2','DPPA5','HESX1','CDX2','GSC','BMP4','ID2','NOG',	'ENSCJAG00000013903','CHRD','ENSCJAG00000014925','ID2',	'ID1','WNT5B','HESX1','DNMT3B','CDX2','TBX3','HAND1','DAZL','MAEL','PRAME','KLF4','PRDM1','BMP4','BMP6','GC','PTCH2','GLI1'})	

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2b,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');
[OBJ3b,a3,b3,OutputCS7b] = transformCS7(OBJ3,OutputCS7,'all');
[OBJ2_2b,a2_2,b2_2,OutputCS62b] = transformCS62(OBJ2_2,OutputCS62,'all');

[OBJ1c,c1,d1,OutputCS5c] = transformCS5(OBJ1,OutputCS5,'notall');
[OBJ2c,c2,d2,OutputCS6c] = transformCS6(OBJ2,OutputCS6,'notall');
[OBJ3c,c3,d3,OutputCS7c] = transformCS7(OBJ3,OutputCS7,'notall');
[OBJ2_2c,c2_2,d2_2,OutputCS62c] = transformCS62(OBJ2_2,OutputCS62,'notall');


%Find where the points are
sh1 = OutputCS5c.cleanShot;
sh2 = OutputCS6c.cleanShot;
sh3 = OutputCS62c.cleanShot;
sh4 = OutputCS7c.cleanShot;

for i = 1:length(sh2)
    sh2{i} = ['E15C2_' sh2{i}];
end

clear ind1
X1 = zeros(length(sh1),4);
for i = 1:length(sh1)
    try
    ind1(i,1) = find(strcmp(sh1{i},Location)==1); 
    X1(i,1:3) = OutputCS5c.cleanX(i,:);
    X1(i,4) = Cl(ind1(i,1));
    catch
    ind1(i,1) = NaN; 
    X1(i,:) = NaN;
    end
end

clear ind1
X2 = zeros(length(sh2),4);
for i = 1:length(sh2)
    try
    ind1(i,1) = find(strcmp(sh2{i},Location)==1); 
    X2(i,1:3) = OutputCS6c.cleanX(i,:);
    X2(i,4) = Cl(ind1(i,1));
    catch
    ind1(i,1) = NaN; 
    X2(i,:) = NaN;
    end
end

clear ind1
X3 = zeros(length(sh3),4);
for i = 1:length(sh3)
    try
    ind1(i,1) = find(strcmp(sh3{i},Location)==1); 
    X3(i,1:3) = OutputCS62c.cleanX(i,:);
    X3(i,4) = Cl(ind1(i,1));
    catch
    ind1(i,1) = NaN; 
    X3(i,:) = NaN;
    end
end

clear ind1
X4 = zeros(length(sh4),4);
for i = 1:length(sh4)
    try
    ind1(i,1) = find(strcmp(sh4{i},Location)==1); 
    X4(i,1:3) = OutputCS7c.cleanX(i,:);
    X4(i,4) = Cl(ind1(i,1));
    catch
    ind1(i,1) = NaN; 
    X4(i,:) = NaN;
    end
end


opac = 0.1;

i = 100

for i = unique(Cl)'


h1=subplot(1,1,1)
indT = find(X1(:,4)==i);
f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices',  OBJ1c.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ1c.objects(20).data.vertices,'Vertices',OBJ1c.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',opac);
%f3 = patch('Faces',OBJ1c.objects(8).data.vertices,'Vertices',OBJ1c.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',opac);
hold on
scatter3(X1(indT,1),X1(indT,2),X1(indT,3),10,'ko','fill')
axis equal
axis off
view([138.5913,5.2765])
camlight('left')
view([c1,d1])
material dull 
colormap(parula);
h1.XLim = h1.XLim*1.2;
h1.YLim = h1.YLim*1.2;

print(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS5_cluster' num2str(i)],'-dpdf','-r1000');
savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS5_cluster' num2str(i) '.fig'])
clf

h2=subplot(1,1,1)
indT = find(X2(:,4)==i);
f1 = patch('Faces',OBJ2c.objects(20).data.vertices,'Vertices',  OBJ2c.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',OBJ2c.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',opac);
f3 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',OBJ2c.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',opac);
f4 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',OBJ2c.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',opac);
hold on
scatter3(X2(indT,1),X2(indT,2),X2(indT,3),10,'ko','fill')
axis equal
axis off
view([8.5335,38.4802])
camlight('left')
view([c2,d2])
%ax = gca(h3)
h2.XLim = h2.XLim*1.2;
h2.YLim = h2.YLim*1.2;
material dull 
colormap(parula);



print(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS6_cluster' num2str(i)],'-dpdf','-r1000');
savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS6_cluster' num2str(i) '.fig'])
clf

h3=subplot(1,1,1)
indT = find(X3(:,4)==i);
f1 = patch('Faces',OBJ2_2c.objects(4).data.vertices,'Vertices',  OBJ2_2c.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ2_2c.objects(8).data.vertices,'Vertices',OBJ2_2c.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',opac);
f4 = patch('Faces',OBJ2_2c.objects(20).data.vertices,'Vertices',OBJ2_2c.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',opac);
hold on
scatter3(X3(indT,1),X3(indT,2),X3(indT,3),10,'ko','fill')
axis equal
axis off
view([c2_2,d2_2])
camlight('left')
material dull 
colormap(parula);


print(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS62_cluster' num2str(i)],'-dpdf','-r1000');
savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/Cl62_cluster' num2str(i) '.fig'])
clf
% 
% h4=subplot(1,1,1)
% indT = find(X4(:,4)==i);
% f1 = patch('Faces',OBJ3c.objects(16).data.vertices,'Vertices', OBJ3c.vertices,'FaceColor',[26,8,115]/255,'LineStyle','none','FaceAlpha',opac);
% f2 = patch('Faces',OBJ3c.objects(20).data.vertices,'Vertices',OBJ3c.vertices,'FaceColor',[2,51,191]/255,'LineStyle','none','FaceAlpha',opac);
% f8 = patch('Faces',OBJ3c.objects(30).data.vertices,'Vertices',OBJ3c.vertices,'FaceColor',[96, 56, 19]/255,'LineStyle','none','FaceAlpha',opac);
% hold on
% scatter3(X4(indT,1),X4(indT,2),X4(indT,3),10,'ko','fill')
% axis equal
% axis off
% view([77.0455, -27.1350])
% camlight('left')
% view([c3,d3])
% material dull 
% colormap(parula);
% 
% 
% print(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS7_cluster' num2str(i)],'-dpdf','-r1000');
% savefig(['/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/NATURE REVISION/Spatial models for revisions/InferredExpressionPatterns/CS7_cluster' num2str(i) '.fig'])
% clf

end
