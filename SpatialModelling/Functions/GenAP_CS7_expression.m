clear all
addpath(genpath('./Functions'))

%Load in the 3D scaffold
[OBJ1,section] = LoadCS7('3D');
%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS7');
[Output] = loadCS7Scaffold(D,Locations,Shots);

load(['~/Desktop/CS7_Hyp1_batch1.mat'],'Hyp1')
load(['~/Desktop/CS7_Hyp2_batch1.mat'],'Hyp2')      
load('./Data/CS7_EmDisc.mat')

[Output] = MarmosetGP_CS7_v3Opt(D,Output,1);


X = zeros(32323,2);
for i = 1:32323
        if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 32323'])
        end

    
[Output] = MarmosetGPInfer_CS7_v3Opt(Output,Line,i,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m2(1);
X(i,2) = Output.m2(end);

end

save(['~/Desktop/CS7_X1_batch1.mat'],'X')
        
return 
X = zeros(10000,2);

load(['~/Desktop/CS7_Hyp1_batch2.mat'],'Hyp1')
load(['~/Desktop/CS7_Hyp2_batch2.mat'],'Hyp2')        
for i = 1:10000
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 10000'])
    end
    
[Output] = MarmosetGPInfer_CS7_v3Opt(Output,Line,i+10000,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m2(1);
X(i,2) = Output.m2(end);
    
end
save(['~/Desktop/CS7_X1_batch2.mat'],'X')

X = zeros(12323,2);

load(['~/Desktop/CS7_Hyp1_batch3.mat'],'Hyp1')
load(['~/Desktop/CS7_Hyp2_batch3.mat'],'Hyp2')   

for i = 1:12323
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 12323'])
    end
    
[Output] = MarmosetGPInfer_CS7_v3Opt(Output,Line,i+20000,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m2(1);
X(i,2) = Output.m2(end);    
 
end
save(['~/Desktop/CS7_X1_batch3.mat'],'X')

              
