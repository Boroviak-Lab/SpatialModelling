function Output = virtualTOMO_seq_CS5(Output)

ind1 = find(contains(Output.cleanShot,'_193_')==1 & strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1 );
ind2 = find(contains(Output.cleanShot,'_198_')==1 & strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1 );
ind3 = find(contains(Output.cleanShot,'_204_')==1 & strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1 );
ind4 = find(contains(Output.cleanShot,'_207_')==1 & strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1 );
ind5 = find(contains(Output.cleanShot,'_213_')==1 & strcmp(Output.cleanAnotaton,'EmDisc_CS5')==1 );


Yts(1) = mean(Output.Ytrain(ind1));
Yts(2) = mean(Output.Ytrain(ind2));
Yts(3) = mean(Output.Ytrain(ind3));
Yts(4) = mean(Output.Ytrain(ind4));
Yts(5) = mean(Output.Ytrain(ind5));

Ytsstd(1) = std(Output.Ytrain(ind1));
Ytsstd(2) = std(Output.Ytrain(ind2));
Ytsstd(3) = std(Output.Ytrain(ind3));
Ytsstd(4) = std(Output.Ytrain(ind4));
Ytsstd(5) = std(Output.Ytrain(ind5));

Output.Yts = Yts;
Output.Ytsstd = Ytsstd;