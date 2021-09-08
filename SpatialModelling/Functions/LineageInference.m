addpath(genpath('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/code/gpml-matlab-v3.6-2015-07-07'))
addpath(genpath('Functions'))
%Load CS7
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,ShotsCS5] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,ShotsCS5);
%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');
[OutputCS6] = loadCS6Scaffold(D,Locations,ShotsCS6);

%Load CS7
%[OBJ3,section] = LoadCS7('3D');
%[D,Locations,XYZ,CellType,ShotsCS7] = LoadShots('CS7');

%These are the subgroups we have
%out =  strrep( ID1.textdata(1,1),'"','');
%out = regexp(out, ',', 'split');
%subpops = out{1};

From = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/LineageInference/Idents2.csv')

S1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/LineageInference/S1_CS5_3k.csv')
C1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/LineageInference/C1_CS5_3k.csv')

LOCS = C1.textdata(2:end,12:21);

shots = ShotsCS6.textdata(2:end,1);
X = NaN*ones(size(LOCS,1),size(LOCS,2));
Y = NaN*ones(size(LOCS,1),size(LOCS,2));
Z = NaN*ones(size(LOCS,1),size(LOCS,2));

for i = 1:length(shots)
    shots2{i} = ['E15C2_' shots{i}];
end

shotssoruce = ShotsCS5.textdata(2:end,1);

%Target locations
shotindex = zeros(length(shots),1);
for i = 1:size(LOCS,1)
    for j = 1:size(LOCS,2)    
    try
       shotindex = find(strcmp(shots2,LOCS{i,j})==1);       
       X(i,j) = ShotsCS6.data(shotindex,1);
       Y(i,j) = ShotsCS6.data(shotindex,2);
       Z(i,j) = ShotsCS6.data(shotindex,3);       
    end
    end
end

%Input locations
Xin = NaN*ones(length(From.textdata(2:end,3)),1);
Yin = NaN*ones(length(From.textdata(2:end,3)),1);
Zin = NaN*ones(length(From.textdata(2:end,3)),1);
for i = 1:length(From.textdata(2:end,3))
    try
        shotindex = find(strcmp(shotssoruce,From.textdata(i+1,3))==1);
       Xin(i,1) = ShotsCS5.data(shotindex,1);
       Yin(i,1) = ShotsCS5.data(shotindex,2);
       Zin(i,1) = ShotsCS5.data(shotindex,3);        
    end    
end

newX = NaN*ones(size(LOCS,1),1);
newY = NaN*ones(size(LOCS,1),1);
newZ = NaN*ones(size(LOCS,1),1);

for i = 1:size(X,1)
    
    u = S1.data(i,:);
    iX = X(i,:);
    
    inds = find(isnan(iX)==0 & u~=0);
    if isempty(inds)==0
        v = u(inds);
        v = v./sum(v);
        newX(i,1) = sum(v.*X(i,inds));
        newY(i,1) =sum(v.*Y(i,inds));
        newZ(i,1) =sum(v.*Z(i,inds));                
    end
    
end

Types = From.textdata(2:end,2);

mapped = find(isnan(newX)==0);
mappedX = [newX(mapped),newY(mapped),newZ(mapped)];
sourceX = [Xin(mapped),Yin(mapped),Zin(mapped)];
Types = Types(mapped);
mapped2 = find(isnan(sourceX(:,1))==0);
mappedX = mappedX(mapped2,:);
sourceX = sourceX(mapped2,:);
Types = Types(mapped2);


OutputCS5.cleanX = sourceX;
OutputCS6.cleanX = mappedX;

[OBJ1_2,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');
[OBJ2_2,a2,b2,OutputCS6b] = transformCS6(OBJ2,OutputCS6,'all');


inds1 = find(strcmp(Types,'Am_CS5'))

for i = inds1'
h = figure(1)
subplot(1,2,1)
opac = .1;
f1 = patch('Faces',OBJ1_2.objects(4).data.vertices,'Vertices',  OBJ1_2.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ1_2.objects(20).data.vertices,'Vertices',OBJ1_2.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',opac);
%f3 = patch('Faces',OBJ1_2.objects(8).data.vertices,'Vertices',OBJ1_2.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',opac);
%f4 = patch('Faces',OBJ1_2.objects(12).data.vertices,'Vertices',OBJ1_2.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',opac);
%f5 = patch('Faces',OBJ1_2.objects(24).data.vertices,'Vertices',OBJ1_2.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',opac);
%f6 = patch('Faces',OBJ1_2.objects(16).data.vertices,'Vertices',OBJ1_2.vertices,'FaceColor',[230,200,0]/255,'LineStyle','none','FaceAlpha',opac);
hold on
scatter3(OutputCS5b.cleanX(i,1),OutputCS5b.cleanX(i,2),OutputCS5b.cleanX(i,3),50,'r','filled')
view([a1,b1])

axis equal
axis off
camlight('right')
material dull 
colormap(parula);


subplot(1,2,2)
f1 = patch('Faces',OBJ2_2.objects(20).data.vertices,'Vertices',  OBJ2_2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ2_2.objects(4).data.vertices,'Vertices',OBJ2_2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',opac);
%f3 = patch('Faces',OBJ2_2.objects(8).data.vertices,'Vertices',OBJ2_2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',opac);
%f4 = patch('Faces',OBJ2_2.objects(24).data.vertices,'Vertices',OBJ2_2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',opac);
%f5 = patch('Faces',OBJ2_2.objects(16).data.vertices,'Vertices',OBJ2_2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',opac);
%f6 = patch('Faces',OBJ2_2.objects(12).data.vertices,'Vertices',OBJ2_2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',opac);
%f7 = patch('Faces',OBJ2_2.objects(28).data.vertices,'Vertices',OBJ2_2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',opac);
%f8 = patch('Faces',OBJ2_2.objects(32).data.vertices,'Vertices',OBJ2_2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',opac);
hold on
scatter3(OutputCS6b.cleanX(i,1),OutputCS6b.cleanX(i,2),OutputCS6b.cleanX(i,3),50,'r','filled')
axis equal
axis off
view([-10.4139,-17.3976])
camlight('left')
material dull 
colormap(parula);

pause
clf
end