function [Output] = MarmosetGPInfer_CS6_v2(Output,OBJ);

%Scale the input data
Xtrain = Output.cleanX / Output.scalefactor;
Ytrain = Output.Ytrain;

HYP2 = Output.HYP1;
HYP3 = Output.HYP2;

im1 = Output.im1;
im2 = Output.im2;

%First is all MODEL 8
par = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain'};
hyp2 = HYP2{1,8};
hyp3 = HYP3{1,8};
[m8 s8] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain', OBJ.vertices/Output.scalefactor);
%[m8_1 s8_1] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain', OBJ.vertices/400);
    
%First is amnion
indT = find(Output.cleanAnotaton=="Am_CS6");
hyp2 = HYP2{1,1};
hyp3 = HYP3{1,1};
[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);

%Second is EmDisk
indT = find(Output.cleanAnotaton=="EmDisc_CS6");
hyp2 = HYP2{1,2};
hyp3 = HYP3{1,2};
[m2 s2] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);

%Second is EmDisk
indT = find(Output.cleanAnotaton=="PGC_CS6");
if isempty(indT)==1
    %Embryo has no PGC
    m3 = [];
else
hyp2 = HYP2{1,3};
hyp3 = HYP3{1,3};
[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);
end

%Third is VE
indT = find(Output.cleanAnotaton=="VE_CS6");
hyp2 = HYP2{1,4};
hyp3 = HYP3{1,4};
[m4 s4] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);


%Next we have SYS
indT = find(Output.cleanAnotaton=="SYS_CS6");
hyp2 = HYP2{1,5};
hyp3 = HYP3{1,5};
[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);

%ExMes
indT = find(Output.cleanAnotaton=="ExMes_CS6" | Output.cleanAnotaton=="Stalk_CS6");
hyp2 = HYP2{1,6};
hyp3 = HYP3{1,6};
[m6 s6] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);
    
%Tb
indT = find(Output.cleanAnotaton=="Tb_CS6");
hyp2 = HYP2{1,7};
hyp3 = HYP3{1,7};
[m7 s7] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/400);
 
Output.m1 = m1;
Output.m2 = m2;
Output.m3 = m3;
Output.m4 = m4;
Output.m5 = m5;
Output.m6 = m6;
Output.m7 = m7;
Output.m8 = m8;

