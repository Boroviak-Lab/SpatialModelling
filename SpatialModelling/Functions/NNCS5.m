function Output = NNCS5(Output,NNs,CynoData,IDcyno,gene)

cynoType = IDcyno.textdata(2:end,2);

mapID = NNs.textdata(2:end,3);
mapID = erase(mapID,"P1_");
mapID = erase(mapID,"P2_");
mapped = NNs.textdata(2:end,4);
score = NNs.data(:,1);

headernames = CynoData.textdata(1,2:end);
headernames = regexp(headernames,'[^""]*','match','once');

genenames = CynoData.textdata(2:end,1);

YtrainCyno = Output.Ytrain*NaN;

for i = 1:length(Output.cleanShot)
    iind = find(strcmp(Output.cleanShot{i}, mapID)==1);
    targmapped = mapped{iind};
    scoremapped = score(iind);
     if scoremapped>0
        cellind = find(strcmp(headernames,targmapped)==1);
        geneind = find(strcmp(genenames,gene)==1);
        YtrainCyno(i) = CynoData.data(geneind,cellind);
        %newAnno{i} = cynoType(iind);
     end
end


Output.fullYcynotrain = YtrainCyno;

%newAnno = newAnno( find(~isnan(YtrainCyno)) );
cynoAnotaton = Output.cleanAnotaton( find(~isnan(YtrainCyno)) );
cynoX = Output.cleanX( find(~isnan(YtrainCyno)) , : );
Ycynotrain = YtrainCyno( find(~isnan(YtrainCyno)) );

%Output.newAnno = newAnno;
Output.cynoAnotaton = cynoAnotaton;
Output.cynoX = cynoX;
Output.Ycynotrain = Ycynotrain;

