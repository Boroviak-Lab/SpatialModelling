addpath(genpath('Functions/'))

D3 = importdata('Data/Keycorrect.csv')
Dexp = importdata('Data/NormData.csv');

[DexpCS5,Arc,Xtest,OutputCS5]  = loadCS5(D3,Dexp);
[DexpCS6,Arc,Xtest,OutputCS6]  = loadCS6(D3,Dexp);
[DexpCS7,Arc,Xtest,OutputCS7]  = loadCS7(D3,Dexp);

Arc1 = importdata('Data/EmDiscCS5_Arc.mat');
Arc2 = importdata('Data/EmDiscCS6_Arc.mat');
Arc3 = importdata('Data/EmDiscCS7_Arc.mat');

Dlist = {'POU5F1','PRDM14','T'}

for j = 1:length(Dlist)
        try
          
[m_CS5,S_CS5,Hyp_CS5] = Marmoset3D_CS5EmDisc_v2(DexpCS5,OutputCS5,'Base2',Arc1,Arc1,Dlist{j});
[m_CS6,S_CS6,Hyp_CS6] = Marmoset3D_CS6EmDisc_v2(DexpCS6,OutputCS6,'Base2',Arc2,Arc1,Dlist{j});
[m_CS7,S_CS7,Hyp_CS7] = Marmoset3D_CS7EmDisc_v2(DexpCS7,OutputCS7,'Base2',Arc3,Dlist{j});

h=figure(j)
subplot(4,6,[1 2 7 8]);
scatter3(OutputCS5.Xtest(OutputCS5.ind1,1),OutputCS5.Xtest(OutputCS5.ind1,2),OutputCS5.Xtest(OutputCS5.ind1,3),25,m_CS5(OutputCS5.ind1),'fill')
hold on
scatter3(OutputCS5.Xtest(OutputCS5.ind4,1),OutputCS5.Xtest(OutputCS5.ind4,2),OutputCS5.Xtest(OutputCS5.ind4,3),25,m_CS5(OutputCS5.ind4),'fill')
view([ -79.1408,77.6211])
set(gca,'clim',[min(min([m_CS5;m_CS6;m_CS7])),max(max([m_CS5;m_CS6;m_CS7]))] )
yl = ylim;
xl = xlim;
zl = zlim;
set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
set(gca,'fontsize', 12)
set(gca,'TickLength',[0 0])
title([Dlist{j}])
set(gca,'linewidth',.7)


subplot(4,6,[3 4 9 10]);
scatter3(OutputCS6.Xtest(OutputCS6.ind1,1),OutputCS6.Xtest(OutputCS6.ind1,2),OutputCS6.Xtest(OutputCS6.ind1,3),25,m_CS6(OutputCS6.ind1),'fill')
hold on
scatter3(OutputCS6.Xtest(OutputCS6.ind4,1),OutputCS6.Xtest(OutputCS6.ind4,2),OutputCS6.Xtest(OutputCS6.ind4,3),25,m_CS6(OutputCS6.ind4),'fill')
view([11.3486,74.1824])
set(gca,'clim',[min(min([m_CS5;m_CS6;m_CS7])),max(max([m_CS5;m_CS6;m_CS7]))] )
yl = ylim;
xl = xlim;
zl = zlim;
set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
set(gca,'TickLength',[0 0])
set(gca,'fontsize', 12)
set(gca,'linewidth',.7)


subplot(4,6,[5 6 11 12]);
scatter3(OutputCS7.Xtest(OutputCS7.ind1,1),OutputCS7.Xtest(OutputCS7.ind1,2),OutputCS7.Xtest(OutputCS7.ind1,3),25,m_CS7(OutputCS7.ind1),'fill')
hold on
scatter3(OutputCS7.Xtest(OutputCS7.ind3,1),OutputCS7.Xtest(OutputCS7.ind3,2),OutputCS7.Xtest(OutputCS7.ind3,3),25,m_CS7(OutputCS7.ind3),'fill')
view([11.3486,74.1824])
set(gca,'clim',[min(min([m_CS5;m_CS6;m_CS7])),max(max([m_CS5;m_CS6;m_CS7]))] )
yl = ylim;
xl = xlim;
zl = zlim;
set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
set(gca,'TickLength',[0 0])
set(gca,'linewidth',.7)


cb=colorbar;
cb.Position = cb.Position + [0.1 -0.02 0 0]
cy=get(cb,'YTick');
set(cb,'YTick',[cy(1),cy(3),cy(end)])
set(gca,'fontsize', 12)

z1 = linspace(0,1,1000)';

m1a = Hyp_CS5.marc1;
s1a = Hyp_CS5.sarc1;
f1a = [m1a+sqrt(s1a); flipdim(m1a-sqrt(s1a),1)];

m2a = Hyp_CS5.marc2;
s2a = Hyp_CS5.sarc2;
f2a = [m2a+sqrt(s2a); flipdim(m2a-sqrt(s2a),1)];

m1b = Hyp_CS6.marc1;
s1b = Hyp_CS6.sarc1;
f1b = [m1b+sqrt(s1b); flipdim(m1b-sqrt(s1b),1)];

m2b = Hyp_CS6.marc2;
s2b = Hyp_CS6.sarc2;
f2b = [m2b+sqrt(s2b); flipdim(m2b-sqrt(s2b),1)];

m1c = Hyp_CS7.marc1;
s1c = Hyp_CS7.sarc1;
f1c = [m1c+sqrt(s1c); flipdim(m1c-sqrt(s1c),1)];

m2c = Hyp_CS7.marc2;
s2c = Hyp_CS7.sarc2;
f2c = [m2c+sqrt(s2c); flipdim(m2c-sqrt(s2c),1)];

subplot(4,6,[13 14]);
fill([z1; flipdim(z1,1)], f1a, [6 6 6]/8);
hold on
plot(z1,m1a,'k','LineWidth',1)
set(gca,'XTick',[])
ylim([min([f1a;f2a]),max([f1a;f2a])])
yt = yticks;
if length(yt)==5,yt = yt(1:2:end),set(gca,'YTick',yt); end
title(['LR: ' num2str(round(-Hyp_CS5.L1(2) -- Hyp_CS5.L2(2),2)) '       '])
set(gca,'fontsize', 8)

subplot(4,6,[19 20]);
fill([z1; flipdim(z1,1)], f2a, [6 6 6]/8);
hold on
plot(z1,m2a,'k','LineWidth',1)
set(gca,'XTick',[])
ylim([min([f1a;f2a]),max([f1a;f2a])])
yt = yticks;
if length(yt)==5,yt = yt(1:2:end),set(gca,'YTick',yt); end
set(gca,'fontsize', 8)


subplot(4,6,[15 16]);
fill([z1; flipdim(z1,1)], f1b, [6 6 6]/8);
hold on
plot(z1,m1b,'k','LineWidth',1)
set(gca,'XTick',[])
ylim([min([f1b;f2b]),max([f1b;f2b])])
yt = yticks;
if length(yt)==5,yt = yt(1:2:end),set(gca,'YTick',yt); end
title(['LR: ' num2str(round(-Hyp_CS6.L1(2) -- Hyp_CS6.L2(2),2)) '       '])
set(gca,'fontsize', 8)

subplot(4,6,[21 22]);
fill([z1; flipdim(z1,1)], f2b, [6 6 6]/8);
hold on
plot(z1,m2b,'k','LineWidth',1)
set(gca,'XTick',[])
ylim([min([f1b;f2b]),max([f1b;f2b])])
yt = yticks;
if length(yt)==5,yt = yt(1:2:end),set(gca,'YTick',yt); end
set(gca,'fontsize', 8)

subplot(4,6,[17 18]);
fill([z1; flipdim(z1,1)], f1c, [6 6 6]/8);
hold on
plot(z1,m1c,'k','LineWidth',1)
set(gca,'XTick',[])
ylim([min([f1c;f2c]),max([f1c;f2c])])
set(gca,'fontsize', 8)
yt = yticks;
if length(yt)==5,yt = yt(1:2:end),set(gca,'YTick',yt); end
title(['LR: ' num2str(round(-Hyp_CS7.L1(3) -- Hyp_CS7.L2(3),2)) '       '])
set(gca,'fontsize', 8)

subplot(4,6,[23 24]);
fill([z1; flipdim(z1,1)], f2c, [6 6 6]/8);
hold on
plot(z1,m2c,'k','LineWidth',1)
set(gca,'XTick',[])
ylim([min([f1c;f2c]),max([f1c;f2c])])
yt = yticks;
if length(yt)==5,yt = yt(1:2:end),set(gca,'YTick',yt); end
set(gca,'fontsize', 8)

%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15.98/2 11.93/2.7])
%ßprint(h, '-dpsc2',['~/Desktop/APAxisCS567_' Dlist{j}],'-r500','-painters')

%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15.98/2 11.93/2.7])
%print(h, '-dpsc2',['~/Desktop/APAxisCS567_' Dlist{j}],'-r500','-painters')

        catch
        end

end