
addpath(genpath('./Functions'))

%Load in the 3D scaffold
[OBJ1,section] = LoadCS5('3D');
%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS5');
[Output] = loadCS5Scaffold(D,Locations,Shots);

L1 = zeros(32323,8);
L2 = zeros(32323,8);
Hyp1 = cell(32323,8);
Hyp2 = cell(32323,8);

load(['~/Desktop/L1_batch1.mat'],'L1')
load(['~/Desktop/L2_batch1.mat'],'L2')
load(['~/Desktop/Hyp1_batch1.mat'],'Hyp1')
load(['~/Desktop/Hyp2_batch1.mat'],'Hyp2')  
        


for i = 6000:32323
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 32323'])
        save(['~/Desktop/L1_batch1.mat'],'L1')
        save(['~/Desktop/L2_batch1.mat'],'L2')
        save(['~/Desktop/Hyp1_batch1.mat'],'Hyp1')
        save(['~/Desktop/Hyp2_batch1.mat'],'Hyp2')        
    end
    
[Output] = MarmosetGP_CS5Opt(D,Output,i);
L1(i,:)= Output.L1;
L2(i,:)= Output.L2;
Hyp1(i,:) = Output.HYP1;
Hyp2(i,:) = Output.HYP2;

end

save(['~/Desktop/L1_batch1.mat'],'L1')
save(['~/Desktop/L2_batch1.mat'],'L2')
save(['~/Desktop/Hyp1_batch1.mat'],'Hyp1')
save(['~/Desktop/Hyp2_batch1.mat'],'Hyp2')  
        
%         
% sdasdsdsada
%         
% 
% 
% L1 = zeros(10000,8);
% L2 = zeros(10000,8);
% Hyp1 = cell(10000,8);
% Hyp2 = cell(10000,8);
% 
% for i = 1:12323
%     if int64(i/1000)==(i/1000)
%         disp(['Step ' num2str(i) ' of 20000'])
%         save(['~/Desktop/L1_batch2.mat'],'L1')
%         save(['~/Desktop/L2_batch2.mat'],'L2')
%         save(['~/Desktop/Hyp1_batch2.mat'],'Hyp1')
%         save(['~/Desktop/Hyp2_batch2.mat'],'Hyp2')        
%     end
%     
% [Output] = MarmosetGP_CS5Opt(D,Output,i+10000);
% L1(i,:)= Output.L1;
% L2(i,:)= Output.L2;
% Hyp1(i,:) = Output.HYP1;
% Hyp2(i,:) = Output.HYP2;
% 
% end
% 
% save(['~/Desktop/L1_batch2.mat'],'L1')
% save(['~/Desktop/L3_batch2.mat'],'L2')
% save(['~/Desktop/Hyp1_batch2.mat'],'Hyp1')
% save(['~/Desktop/Hyp2_batch2.mat'],'Hyp2')  
%         
%                 
% 
% L1 = zeros(12323,8);
% L2 = zeros(12323,8);
% Hyp1 = cell(12323,8);
% Hyp2 = cell(12323,8);
% 
% for i = 1:12323
%     if int64(i/1000)==(i/1000)
%         disp(['Step ' num2str(i) ' of 20000'])
%         save(['~/Desktop/L1_batch3.mat'],'L1')
%         save(['~/Desktop/L2_batch3.mat'],'L2')
%         save(['~/Desktop/Hyp1_batch3.mat'],'Hyp1')
%         save(['~/Desktop/Hyp2_batch3.mat'],'Hyp2')        
%     end
%     
% [Output] = MarmosetGP_CS5Opt(D,Output,i+20000);
% L1(i,:)= Output.L1;
% L2(i,:)= Output.L2;
% Hyp1(i,:) = Output.HYP1;
% Hyp2(i,:) = Output.HYP2;
% 
% end
% 
% save(['~/Desktop/L1_batch3.mat'],'L1')
% save(['~/Desktop/L3_batch3.mat'],'L2')
% save(['~/Desktop/Hyp1_batch3.mat'],'Hyp1')
% save(['~/Desktop/Hyp2_batch3.mat'],'Hyp2')  
%         
% 
% 
% 
% 
% 
% L1 = zeros(12323,8);
% L2 = zeros(12323,8);
% Hyp1 = cell(12323,8);
% Hyp2 = cell(12323,8);
% 
% for i = 1:12323
%     if int64(i/1000)==(i/1000)
%         disp(['Step ' num2str(i) ' of 20000'])
%         save(['~/Desktop/L1_batch3.mat'],'L1')
%         save(['~/Desktop/L2_batch3.mat'],'L2')
%         save(['~/Desktop/Hyp1_batch3.mat'],'Hyp1')
%         save(['~/Desktop/Hyp2_batch3.mat'],'Hyp2')        
%     end
%     
% [Output] = MarmosetGP_CS5Opt(D,Output,i+20000);
% L1(i,:)= Output.L1;
% L2(i,:)= Output.L2;
% Hyp1(i,:) = Output.HYP1;
% Hyp2(i,:) = Output.HYP2;
% 
% end
% 
% save(['~/Desktop/L1_batch4.mat'],'L1')
% save(['~/Desktop/L3_batch4.mat'],'L2')
% save(['~/Desktop/Hyp1_batch4.mat'],'Hyp1')
% save(['~/Desktop/Hyp2_batch4.mat'],'Hyp2')  
%         
%               
              
