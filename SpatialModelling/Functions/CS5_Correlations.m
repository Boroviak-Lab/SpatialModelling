addpath(genpath('./Functions'))

%Load in the 3D scaffold
[OBJ1,section] = LoadCS5('3D');

%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');

%Process the shots onto scaffold
[Output] = loadCS5Scaffold(D,Locations,Shots);

Ytrain = D.data(:,Output.cleanindex);
size(Output.Xtrain)

C = corr(Ytrain);
d = dist(Output.cleanX');
U = triu(C,1);
d = triu(d,1);
U(find(d==0)) = NaN;
d(find(d==0)) = NaN;

C1=(d(~isnan(U)));
U1=(U(~isnan(U)));

subplot(2,2,1);
plot(C1,U1,'.')

subplot(2,2,2);
ind = find(strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1)
U1 = U(ind,ind);
d1 = d(ind,ind);
C1=(d1(~isnan(U1)));
U1=(U1(~isnan(U1)));
plot(C1,U1,'.')

subplot(2,2,3);
ind = find(strcmp(Output.cleanAnotaton,'VE_CS5')==1)
U1 = U(ind,ind);
d1 = d(ind,ind);
C1=(d1(~isnan(U1)));
U1=(U1(~isnan(U1)));
plot(C1,U1,'.')

subplot(2,2,4);
ind = find(strcmp(Output.cleanAnotaton,'Tb_CS5')==1)
U1 = U(ind,ind);
d1 = d(ind,ind);
C1=(d1(~isnan(U1)));
U1=(U1(~isnan(U1)));
plot(C1,U1,'.')




%Seperate by tissue
