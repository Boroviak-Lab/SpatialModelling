addpath(genpath('./Functions'))

%Load in the 3D scaffold
[OBJ1,section] = LoadCS6('3D');
%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS6');
[Output] = loadCS6Scaffold(D,Locations,Shots);


return
%L1 = zeros(32323,9);
%L2 = zeros(32323,9);
%Hyp1 = cell(32323,9);
%Hyp2 = cell(32323,9);

load(['~/Desktop/CS6_L1_batch1.mat'])
load(['~/Desktop/CS6_L2_batch1.mat'])
load(['~/Desktop/CS6_Hyp1_batch1.mat'])
load(['~/Desktop/CS6_Hyp2_batch1.mat'])  

for i = 26000:32323
    if int64(i/1000)==(i/1000)
        disp(['Step ' num2str(i) ' of 32323'])
        save(['~/Desktop/CS6_L1_batch1.mat'],'L1')
        save(['~/Desktop/CS6_L2_batch1.mat'],'L2')
        save(['~/Desktop/CS6_Hyp1_batch1.mat'],'Hyp1')
        save(['~/Desktop/CS6_Hyp2_batch1.mat'],'Hyp2')        
    end
    
[Output] = MarmosetGP_CS6_v3Opt(D,Output,i);
L1(i,:)= Output.L1;
L2(i,:)= Output.L2;
Hyp1(i,:) = Output.HYP1;
Hyp2(i,:) = Output.HYP2;

end

save(['~/Desktop/CS6_L1_batch1.mat'],'L1')
save(['~/Desktop/CS6_L2_batch1.mat'],'L2')
save(['~/Desktop/CS6_Hyp1_batch1.mat'],'Hyp1')
save(['~/Desktop/CS6_Hyp2_batch1.mat'],'Hyp2')  
%         
% L1 = zeros(10000,9);
% L2 = zeros(10000,9);
% Hyp1 = cell(10000,9);
% Hyp2 = cell(10000,9);
% 
% for i = 1:12323
%     if int64(i/1000)==(i/1000)
%         disp(['Step ' num2str(i) ' of 20000'])
%         save(['~/Desktop/CS6_L1_batch2.mat'],'L1')
%         save(['~/Desktop/CS6_L2_batch2.mat'],'L2')
%         save(['~/Desktop/CS6_Hyp1_batch2.mat'],'Hyp1')
%         save(['~/Desktop/CS6_Hyp2_batch2.mat'],'Hyp2')        
%     end
%     
% [Output] = MarmosetGP_CS6_v3Opt(D,Output,i+10000);
% L1(i,:)= Output.L1;
% L2(i,:)= Output.L2;
% Hyp1(i,:) = Output.HYP1;
% Hyp2(i,:) = Output.HYP2;
% 
% end
% 
% save(['~/Desktop/CS6_L1_batch2.mat'],'L1')
% save(['~/Desktop/CS6_L3_batch2.mat'],'L2')
% save(['~/Desktop/CS6_Hyp1_batch2.mat'],'Hyp1')
% save(['~/Desktop/CS6_Hyp2_batch2.mat'],'Hyp2')  
%                       
% L1 = zeros(12323,9);
% L2 = zeros(12323,9);
% Hyp1 = cell(12323,9);
% Hyp2 = cell(12323,9);
% 
% for i = 1:12323
%     if int64(i/1000)==(i/1000)
%         disp(['Step ' num2str(i) ' of 20000'])
%         save(['~/Desktop/CS6_L1_batch3.mat'],'L1')
%         save(['~/Desktop/CS6_L2_batch3.mat'],'L2')
%         save(['~/Desktop/CS6_Hyp1_batch3.mat'],'Hyp1')
%         save(['~/Desktop/CS6_Hyp2_batch3.mat'],'Hyp2')        
%     end
%     
% [Output] = MarmosetGP_CS6_v3Opt(D,Output,i+20000);
% L1(i,:)= Output.L1;
% L2(i,:)= Output.L2;
% Hyp1(i,:) = Output.HYP1;
% Hyp2(i,:) = Output.HYP2;
% 
% end
% 
% save(['~/Desktop/CS6_L1_batch3.mat'],'L1')
% save(['~/Desktop/CS6_L3_batch3.mat'],'L2')
% save(['~/Desktop/CS6_Hyp1_batch3.mat'],'Hyp1')
% save(['~/Desktop/CS6_Hyp2_batch3.mat'],'Hyp2')  
%         
%               
