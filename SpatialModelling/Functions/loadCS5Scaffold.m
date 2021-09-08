function [Output] = loadCS5Scaffold(D,Locations,Shots)

shots = Shots.textdata(2:end,1);

ReferenceLocatioins = Locations.textdata(2:end,2);
ReferenceLocatioins = erase(ReferenceLocatioins,"P1_");
ReferenceLocatioins = erase(ReferenceLocatioins,"P2_");

ReferenceTypes = Locations.textdata(2:end,4);

shotindex = zeros(length(shots),1);
for i = 1:length(shots)
    try
       shotindex(i,1) = find(strcmp(shots{i},ReferenceLocatioins)==1);
    end
end

%Remove the shots that weren't there due to QC issues
QCshotindex = shotindex(find(shotindex~=0));
QCshots = shots(find(shotindex~=0));
XQC = Shots.data(find(shotindex~=0),:);
QCReferenceTypes = ReferenceTypes(QCshotindex);

%Remove mixed shots
Output.mixedindex = QCshotindex;
Output.mixedX = XQC;
Output.mixedShot = QCshots;
Output.mixedAnotaton = QCReferenceTypes;

idex = find( strcmp(QCReferenceTypes,'Am_CS5')==1 | strcmp(QCReferenceTypes,'EmDisc_CS5')==1 | strcmp(QCReferenceTypes,'Am_CS5_PGC')==1 | strcmp(QCReferenceTypes,'EmDisc_CS5')==1 | strcmp(QCReferenceTypes,'ExMes_CS5')==1 | strcmp(QCReferenceTypes,'SYS_CS5')==1 | strcmp(QCReferenceTypes,'Tb_CS5')==1 | strcmp(QCReferenceTypes,'VE_CS5')==1 );

cleanindex = QCshotindex(idex);
cleanshots = QCshots(idex);
cleanX = XQC(idex,:);
cleanReferenceTypes = QCReferenceTypes(idex);

Output.cleanindex = cleanindex;
Output.cleanX = cleanX;
Output.cleanShot = cleanshots;
Output.cleanAnotaton = cleanReferenceTypes;

Output.scalefactor = 400;

%h = scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,'fill','MarkerEdgeColor',[12,156,245]/255, 'MarkerFaceColor',[12,156,245]/255)
%hold on
%scatter3(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3),95,'fill')
%text(Output.Xtrain(find(Output.LabBin==1),1),Output.Xtrain(find(Output.LabBin==1),2),Output.Xtrain(find(Output.LabBin==1),3), CellIDs4(find(Output.LabBin==1)));
%scatter3(Output.Xtest(Output.ind0,1),Output.Xtest(Output.ind0,2),Output.Xtest(Output.ind0,3),5,'fill')
%scatter3(Output.Xtrain(setdiff(find(Output.LabBin==1),indss),1),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),2),Output.Xtrain(setdiff(find(Output.LabBin==1),indss),3),20,'fill','MarkerEdgeColor','k', 'MarkerFaceColor','k')
%text(Output.Xtrain(find(Output.LabBin==0),1),Output.Xtrain(find(Output.LabBin==0),2),Output.Xtrain(find(Output.LabBin==0),3), CellIDs4(find(Output.LabBin==0)));
%scatter3(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3),195,'fill','MarkerEdgeColor','r', 'MarkerFaceColor','r')
%text(Output.Xtrain(indss,1),Output.Xtrain(indss,2),Output.Xtrain(indss,3), CellIDs4(indss));
%view([ -79.1408,77.6211])
%view([-78.2671,26.9624])
%set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
%alpha = 0.05;
%set(h,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
%set(h2,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
