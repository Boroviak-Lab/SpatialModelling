addpath(genpath('./Functions'))

%Load in the 3D scaffold
[OBJ1,section] = LoadCS5('3D');
%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
[Output] = loadCS5Scaffold(D,Locations,Shots);

load(['~/Desktop/Hyp1_batch1.mat'],'Hyp1')
load(['~/Desktop/Hyp2_batch1.mat'],'Hyp2')      

load('./Data/CS5_EmDisc.mat')

[Output] = MarmosetGP_CS5Opt(D,Output,1);

X = zeros(32323,2);

for i = 1:32323
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 32323'])
    end
    

[Output] = MarmosetGPInfer_CS5Opt(Output,Line,i,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m2(1);
X(i,2) = Output.m2(end);    
     
end

save(['~/Desktop/CS5_X1_batch1.mat'],'X')
        

return

load(['~/Desktop/Hyp1_batch2.mat'],'Hyp1')
load(['~/Desktop/Hyp2_batch2.mat'],'Hyp2')      

        
X = zeros(10000,2);

for i = 1:10000
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 20000'])
    end
    

[Output] = MarmosetGPInfer_CS5Opt(Output,Line,i+10000,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m2(1);
X(i,2) = Output.m2(end);    
     
      
    
end

save(['~/Desktop/CS5_X1_batch2.mat'],'X')
        
load(['~/Desktop/Hyp1_batch3.mat'],'Hyp1')
load(['~/Desktop/Hyp2_batch3.mat'],'Hyp2')      
                

X = zeros(12323,2);

for i = 1:12323
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 20000'])
    end
     
[Output] = MarmosetGPInfer_CS5Opt(Output,Line,i+20000,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m2(1);
X(i,2) = Output.m2(end);    
     
  
end

save(['~/Desktop/CS5_X1_batch3.mat'],'X')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%And VE
load(['~/Desktop/Hyp1_batch1.mat'],'Hyp1')
load(['~/Desktop/Hyp2_batch1.mat'],'Hyp2')      

load('./Data/CS5_VE.mat')

X = zeros(10000,2);

for i = 1:10000
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 20000'])
    end
    

[Output] = MarmosetGPInfer_CS5Opt(Output,Line,i,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m4(1);
X(i,2) = Output.m4(end);    
     
end

save(['~/Desktop/CS5_VE_X1_batch1.mat'],'X')
        

load(['~/Desktop/Hyp1_batch2.mat'],'Hyp1')
load(['~/Desktop/Hyp2_batch2.mat'],'Hyp2')      

        
X = zeros(10000,2);

for i = 1:10000
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 20000'])
    end
    

[Output] = MarmosetGPInfer_CS5Opt(Output,Line,i+10000,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m4(1);
X(i,2) = Output.m4(end);    
     
      
    
end

save(['~/Desktop/CS5_VE_X1_batch2.mat'],'X')
        
load(['~/Desktop/Hyp1_batch3.mat'],'Hyp1')
load(['~/Desktop/Hyp2_batch3.mat'],'Hyp2')      
                

X = zeros(12323,2);

for i = 1:12323
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 20000'])
    end
     
[Output] = MarmosetGPInfer_CS5Opt(Output,Line,i+20000,D,Hyp1(i,:),Hyp2(i,:),'Line');

X(i,1) = Output.m4(1);
X(i,2) = Output.m4(end);    
     
end

save(['~/Desktop/CS5_VE_X1_batch3.mat'],'X')
        
              
              
