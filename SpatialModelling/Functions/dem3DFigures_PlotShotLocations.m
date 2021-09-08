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


D1 = importdata('/Users/christopherpenfold/Desktop/Loc.csv')
LOCS = D1.textdata(2:355,2);
IDs  = D1.textdata(2:355,3);
uIDs = unique(unique(D1.textdata(2:355,3)));


pgcshot = intersect(strcat('E15C2_',OutputCS6c.cleanShot(find(strcmp(OutputCS6c.cleanAnotaton,'PGC_CS6')==1))),LOCS)
for i = 1:length(pgcshot)
    try
  cs(i,1) = find(strcmp(pgcshot{i},LOCS)==1);
   catch
  cs(i,1) = 0;
    end
end

for j = 1:length(uIDs)
ind  = find(strcmp(IDs,uIDs{j})==1);

currLOC = LOCS(ind);

for i = 1:length(LOCS)
    try
  cs5ind(i,1) = find(strcmp(currLOC{i},OutputCS5c.cleanShot)==1);
   catch
  cs5ind(i,1) = 0;
    end
end

for i = 1:length(LOCS)
    try
  cs6ind(i,1) = find(strcmp(currLOC{i},strcat('E15C2_',OutputCS6c.cleanShot))==1);
   catch
  cs6ind(i,1) = 0;
    end
end

for i = 1:length(LOCS)
    try
  cs62ind(i,1) = find(strcmp(currLOC{i},OutputCS62c.cleanShot)==1);
   catch
        cs62ind(i,1) = 0;
    end
end

for i = 1:length(LOCS)
    try
  cs7ind(i,1) = find(strcmp(currLOC{i},OutputCS7c.cleanShot)==1);
    catch
        cs7ind(i,1) = 0;
    end
end



%Plot CS5
opac = 0.1;

h0 = figure(1)
h1 = subplot(1,4,1);
f1 = patch('Faces',OBJ1c.objects(4).data.vertices,'Vertices', OBJ1c.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ1c.objects(20).data.vertices,'Vertices',OBJ1c.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',opac);
f3 = patch('Faces',OBJ1c.objects(8).data.vertices,'Vertices', OBJ1c.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',opac);
%f4 = patch('Faces',OBJ1c.objects(12).data.vertices,'Vertices',OBJ1c.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',opac);
%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',opac);
%f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,200,0]/255,'LineStyle','none','FaceAlpha',opac);
axis equal
axis off
camlight('right')
material dull 
colormap(parula);
hold on
scatter3(OutputCS5c.cleanX(cs5ind(find(cs5ind~=0)),1),OutputCS5c.cleanX(cs5ind(find(cs5ind~=0)),2),OutputCS5c.cleanX(cs5ind(find(cs5ind~=0)),3),20,'ko','fill')
view(c1,d1)
h1.XLim = h1.XLim*1.2;
h1.YLim = h1.YLim*1.2;


h2 = subplot(1,4,2);
f1 = patch('Faces',OBJ2c.objects(20).data.vertices,'Vertices',  OBJ2c.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ2c.objects(4).data.vertices,'Vertices',OBJ2c.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',opac);
f3 = patch('Faces',OBJ2c.objects(8).data.vertices,'Vertices',OBJ2c.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',opac);
%f4 = patch('Faces',OBJ2c.objects(24).data.vertices,'Vertices',OBJ2c.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',opac);
%f5 = patch('Faces',OBJ2A.objects(16).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',opac);
%f6 = patch('Faces',OBJ2A.objects(12).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',opac);
%f7 = patch('Faces',OBJ2A.objects(28).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',opac);
%f8 = patch('Faces',OBJ2A.objects(32).data.vertices,'Vertices',OBJ2A.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',opac);
%axis equal
axis off
camlight('right')
material dull 
colormap(parula);
hold on
scatter3(OutputCS6c.cleanX(cs6ind(find(cs6ind~=0)),1),OutputCS6c.cleanX(cs6ind(find(cs6ind~=0)),2),OutputCS6c.cleanX(cs6ind(find(cs6ind~=0)),3),20,'ko','fill')
view(c2,d2)
axis equal
%ax = gca(h4)
h2.XLim = h2.XLim*1.1;
h2.YLim = h2.YLim*1.1;

%h6 = figure(1)
h3 = subplot(1,4,3);
f1 = patch('Faces',OBJ2_2c.objects(4).data.vertices,'Vertices', OBJ2_2c.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ2_2c.objects(8).data.vertices,'Vertices', OBJ2_2c.vertices,'FaceColor',[218,122,7]/255,'LineStyle','none','FaceAlpha',opac);
f4 = patch('Faces',OBJ2_2c.objects(20).data.vertices,'Vertices',OBJ2_2c.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',opac);
%f5 = patch('Faces',OBJ2_2b.objects(12).data.vertices,'Vertices',OBJ2_2b.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',1);
%f6 = patch('Faces',OBJ2_2b.objects(16).data.vertices,'Vertices',OBJ2_2b.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',1);
%f7 = patch('Faces',OBJ2_2b.objects(24).data.vertices,'Vertices',OBJ2_2b.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',1);
%f8 = patch('Faces',OBJ2_2b.objects(28).data.vertices,'Vertices',OBJ2_2b.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',1);
axis equal
axis off
camlight('right')
material dull 
colormap(parula);
hold on
scatter3(OutputCS62c.cleanX(cs62ind(find(cs62ind~=0)),1),OutputCS62c.cleanX(cs62ind(find(cs62ind~=0)),2),OutputCS62c.cleanX(cs62ind(find(cs62ind~=0)),3),20,'ko','fill')
view([c2_2,d2_2])

%h4b = figure(1)
h4 = subplot(1,4,4);
f1 = patch('Faces',OBJ3c.objects(16).data.vertices,'Vertices', OBJ3c.vertices,'FaceColor',[26,8,115]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ3c.objects(20).data.vertices,'Vertices',OBJ3c.vertices,'FaceColor',[2,51,191]/255,'LineStyle','none','FaceAlpha',opac);
f8 = patch('Faces',OBJ3c.objects(30).data.vertices,'Vertices',OBJ3c.vertices,'FaceColor',[96, 56, 19]/255,'LineStyle','none','FaceAlpha',opac);
%f4 = patch('Faces',OBJ3b.objects(12).data.vertices,'Vertices',OBJ3b.vertices,'FaceColor',[191,60,4]/255,'LineStyle','none','FaceAlpha',opac);
%f5 = patch('Faces',OBJ3b.objects(8).data.vertices,'Vertices',OBJ3b.vertices,'FaceColor',[191,113,4]/255,'LineStyle','none','FaceAlpha',opac);
%f6 = patch('Faces',OBJ3b.objects(4).data.vertices,'Vertices',OBJ3b.vertices,'FaceColor',[113,8,166]/255,'LineStyle','none','FaceAlpha',opac);
%f7 = patch('Faces',OBJ3b.objects(34).data.vertices,'Vertices',OBJ3b.vertices,'FaceColor',[150,119,0]/255,'LineStyle','none','FaceAlpha',opac);
axis equal
axis off
camlight('right')
material dull 
colormap(parula);
hold on
scatter3(OutputCS7c.cleanX(cs7ind(find(cs7ind~=0)),1),OutputCS7c.cleanX(cs7ind(find(cs7ind~=0)),2),OutputCS7c.cleanX(cs7ind(find(cs7ind~=0)),3),20,'ko','fill')
view(c3,d3)
h4.XLim = h4.XLim*1.2;
h4.YLim = h4.YLim*1.2;



set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12*0.3942])
print('-dpng',['~/Desktop/ClLocations_' uIDs{j} '.png'],'-r1000')
print('-dpng',['~/Desktop/ClLocations_' uIDs{j} '.pdf'],'-r1000')
%print(h5,['Plots/WholeEmb_' gene{i} '_CS5-7_scaled600.pdf'],'-dpdf','-r1000');

close all

end