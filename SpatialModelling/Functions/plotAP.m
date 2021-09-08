function h2 = plotAP(Output, figno, subplotp, Type);


 h = figure(figno);
 h2 = subplot(subplotp(1),subplotp(2),[subplotp(3:end)]);
 
 if strcmp(Type,'VE')==1
 m = Output.m4;
 s = Output.s4;
 else
 m = Output.m2;
 s = Output.s2;
 end
 z = linspace(0,1,length(m))'; 
 f = [m+sqrt(s); flipdim(m-sqrt(s),1)];
 fill([z; flipdim(z,1)], f, [6 6 6]/8);
hold on 
plot(z,m,'k-','LineWidth',2)
 set(gca,'XTick',[])