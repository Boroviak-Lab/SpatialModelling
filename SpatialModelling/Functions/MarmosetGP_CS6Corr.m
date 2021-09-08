function [Output] = MarmosetGP_CS6Corr(Output);

%Scale the input data
Xtrain = Output.cleanX / Output.scalefactor;

pg1 = {@priorGauss,-0.6931,2};  
pg2 = {@priorGauss,1,2};  
pg3 = {@priorGauss,log(0.1),1};  
pc = {@priorClamped};

%Prior for base model
prior1.mean = {[]};  
prior1.cov  = {pg1;pg2};
prior1.lik  = {pg3};

%prior1.mean = {[]};  
%prior1.cov  = {pg2;pg1};
%prior1.lik  = {pg3};

%Prior for long lenghth-scale model
prior2.mean = {[]};  
prior2.cov  = {pc;pg1};
prior2.lik  = {pg3};

im1 = {@infPrior,@infExact,prior1};
im2 = {@infPrior,@infExact,prior2};
 
HYP2 = cell(size(Output.P,2),8);
HYP3 = cell(size(Output.P,2),8);

L1 = zeros(size(Output.P,2),8);
L2 = zeros(size(Output.P,2),8);

for i = 1:size(Output.P,2)
   
Ytrain = Output.P(:,i)';    
%Initalise a high noiose std and long length scale

hyp1.cov  = [ log(0.1) ; log(std(Ytrain)) ];
hyp1.lik  = [log(std(Ytrain)/4)];    
hyp1.mean = mean(Ytrain);       

hyp1B.cov  = [40 ; log(0.1) ];
hyp1B.lik  = [var(Ytrain)/2];    
hyp1B.mean = mean(Ytrain); 

%First is all
par = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,9} = hyp2;
HYP3{i,9} = hyp3;
%[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,9) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
[L2(i,9) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
    
%First is amnion
indT = find(strcmp(Output.cleanAnotaton,"Am_CS6")==1 | strcmp(Output.cleanAnotaton,"Am_CS6_PGC")==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,1} = hyp2;
HYP3{i,1} = hyp3;
%[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,1) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,1) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


%Second is EmDisk
indT = find(strcmp(Output.cleanAnotaton,"EmDisc_CS6")==1 |  strcmp(Output.cleanAnotaton,"EmDisc_CS6_PGC")==1 | strcmp(Output.cleanAnotaton,"EmDisc_Gast_CS6")==1 | strcmp(Output.cleanAnotaton,"EmDisc_gast_CS6")==1 );
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,2} = hyp2;
HYP3{i,2} = hyp3;
%[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
%[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,2) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,2) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

%Second is PGC
indT = find(strcmp(Output.cleanAnotaton,"PGC_CS6")==1);
    %Embryo 2 has no PGCs
if isempty(indT) == 1
    L1(i,3) = NaN;
    L2(i,3) = NaN;
else
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,3} = hyp2;
HYP3{i,3} = hyp3;
%[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
%[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,3) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,3) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
end

%Third is VE
indT = find(strcmp(Output.cleanAnotaton,"VE_CS6")==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,4} = hyp2;
HYP3{i,4} = hyp3;
%[m2 s2] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,4) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,4) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


%Next we have SYS
indT = find(strcmp(Output.cleanAnotaton,"SYS_CS6")==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,5} = hyp2;
HYP3{i,5} = hyp3;
%[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,5) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,5) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

%ExMes
indT = find(strcmp(Output.cleanAnotaton,"ExMes_CS6")==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,6} = hyp2;
HYP3{i,6} = hyp3;
%[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,6) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,6) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
    
%Tb
indT = find(strcmp(Output.cleanAnotaton,"Tb_CS6")==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,7} = hyp2;
HYP3{i,7} = hyp3;
%[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,7) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,7) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
 


%Second is EmDisk
indT = find(strcmp(Output.cleanAnotaton,"Stalk_CS6")==1 | strcmp(Output.cleanAnotaton,"EmDisc_Stalk_CS6")==1 | strcmp(Output.cleanAnotaton,"EmDisc_stalk_CS6")==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -10000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -10000, im2, par{:});         % optimise
HYP2{i,8} = hyp2;
HYP3{i,8} = hyp3;
%[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
%[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1(i,8) ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2(i,8) ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


end
%Output.m_0 = m0;
%Output.m_1 = m1;
%Output.m_2 = m2;
%Output.m_3 = m3;
%Output.m_5 = m5;



%Output.Ytrain = Ytrain;
%Output.Xtrain = Xtrain;
Output.L1 = L1; %[L1_1,L1_2,L1_3,L1_4,L1_5,L1_6,L1_7,L1_8];
Output.L2 = L2; %[L2_1,L2_2,L2_3,L2_4,L2_5,L2_6,L2_7,L1_8];
Output.HYP1 = HYP2;
Output.HYP2 = HYP3;

Output.im1 = im1;
Output.im2 = im2;
Output.im1 = prior1;
Output.im2 = prior2;

Output.Ytrain = Ytrain;

%Output.hyp2 = hyp1;
%Output.hyp2 = hyp2;
%Output.hyp3 = hyp3;