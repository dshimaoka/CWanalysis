contrastList = [0 .12 .25 .5];
allContrast = (100*[-contrastList contrastList]);

totnEvents = squeeze(nansum(nEvents,3));
ratioRight = squeeze(totnEvents(:,3))./squeeze(sum(totnEvents,2));
ratioLeft = squeeze(totnEvents(:,2))./squeeze(sum(totnEvents,2));
ratioNoGo = squeeze(totnEvents(:,4))./squeeze(sum(totnEvents,2));

tgtConstIdx = [4 3 2 6 7 8]; %1:8;
plot(allContrast(tgtConstIdx), ratioRight(tgtConstIdx),'r');
hold on
plot(allContrast(tgtConstIdx), ratioLeft(tgtConstIdx),'b');
plot(allContrast(tgtConstIdx), ratioNoGo(tgtConstIdx),'g');
legend('choose R','choose L','No Go');
screen2png('psychometricCurve');
