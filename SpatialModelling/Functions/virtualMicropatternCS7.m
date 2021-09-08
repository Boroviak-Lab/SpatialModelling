function [output] = virtualMicropatternCS7(D,Output,Line,genelist,figno)

figure(figno)
%Process the shots onto scaffold
[Output] = MarmosetGP_CS7_v3(D,Output,genelist{1});

%Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS7_v3(Output,Line,'Line');

subplot(2,2,1);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title(genelist{1})
colorbar
axis off

[Output] = MarmosetGP_CS7_v3(D,Output,genelist{2});
%Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS7_v3(Output,Line,'Line');
subplot(2,2,2);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title(genelist{2})
colorbar
axis off

[Output] = MarmosetGP_CS7_v3(D,Output,genelist{3});
%Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS7_v3(Output,Line,'Line');
subplot(2,2,3);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title(genelist{3})
colorbar
axis off

[Output] = MarmosetGP_CS7_v3(D,Output,genelist{4});
%Line.vertices = [x_ad,y_ad,z_ad];
[Output] = MarmosetGPInfer_CS7_v3(Output,Line,'Line');
subplot(2,2,4);
N = 10000;
u = rand(1,N)*2*pi;
ind = randi(length(Output.m2),[1 N] );        
r = normrnd(Output.m2(ind),Output.s2(ind));
scatter(ind.*cos(u),ind.*sin(u),20,r,'filled')
title(genelist{4})
colorbar
axis off
