clear all; close all;

bp     = 7;
mod    = 2;
icend  = 7;
    
plotCatalog(bp,mod,icend);
plotGlobaldat(bp,mod,icend);
plotOffFaultSt(bp,mod,icend);
plotOnFaultSt(bp,mod,icend);
plotSlipStressProfileBP7(bp,mod,icend,1); % mode 1, strike.
plotSlipStressProfileBP7(bp,mod,icend,2); % mode 2, depth.
plotRPT(bp,mod);