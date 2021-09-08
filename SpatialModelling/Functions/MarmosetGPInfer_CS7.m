function [Output] = MarmosetGPInfer_CS7(varargin);


if length(varargin)==2
    Output = varargin{1};
    OBJ = varargin{2};
    savevar = 0;
elseif length(varargin)==3
    Output = varargin{1};
    OBJ = varargin{2};
    savevar = 1;
else
    disp('Wrong number of parameters')
end

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
indT = find(Output.cleanAnotaton=="Am_CS7");
hyp2 = HYP2{1,1};
hyp3 = HYP3{1,1};
[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);

%Second is EmDisk
indT = find(Output.cleanAnotaton=="EmDisc_CS7" | Output.cleanAnotaton=="Stalk_CS7" | Output.cleanAnotaton=="EmDisc_stalk_CS7");
hyp2 = HYP2{1,2};
hyp3 = HYP3{1,2};
[m2 s2] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);

m3 = [];
m4 = [];

%Next we have SYS
indT = find(Output.cleanAnotaton=="SYS_CS7");
hyp2 = HYP2{1,5};
hyp3 = HYP3{1,5};
[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);

%ExMes
indT = find(Output.cleanAnotaton=="ExMes_CS7");
hyp2 = HYP2{1,6};
hyp3 = HYP3{1,6};
[m6 s6] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);
    
%Tb
indT = find(Output.cleanAnotaton=="Tb_CS7");
hyp2 = HYP2{1,7};
hyp3 = HYP3{1,7};
[m7 s7] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ.vertices/Output.scalefactor);
 
Output.m1 = m1;
Output.m2 = m2;
Output.m3 = m3;
Output.m4 = m4;
Output.m5 = m5;
Output.m6 = m6;
Output.m7 = m7;
Output.m8 = m8;

try
Output.cLim1 = [min(Output.m1(OBJ.objects(16).data.vertices)),max(Output.m1(OBJ.objects(16).data.vertices))];
Output.cLim2 = [min(Output.m2(OBJ.objects(20).data.vertices)),max(Output.m2(OBJ.objects(20).data.vertices))];
Output.cLim2b = [min(Output.m2(OBJ.objects(30).data.vertices)),max(Output.m2(OBJ.objects(30).data.vertices))];
Output.cLim5 = [min(Output.m5(OBJ.objects(8).data.vertices)),max(Output.m5(OBJ.objects(8).data.vertices))];
Output.cLim7 = [min(Output.m7(OBJ.objects(4).data.vertices)),max(Output.m7(OBJ.objects(4).data.vertices))];
Output.cLim6 = [min(Output.m6(OBJ.objects(34).data.vertices)),max(Output.m6(OBJ.objects(34).data.vertices))];
Output.cLim = [min( [Output.cLim1,Output.cLim2,Output.cLim2b,Output.cLim5,Output.cLim6,Output.cLim7]), max( [Output.cLim1,Output.cLim2,Output.cLim2b,Output.cLim5,Output.cLim6,Output.cLim7])];
end

if savevar==1
Output.s1 = s1;
Output.s2 = s2;
Output.s3 = [];
Output.s4 = [];
Output.s5 = s5;
Output.s6 = s6;
Output.s7 = s7;
Output.s8 = s8;
end
