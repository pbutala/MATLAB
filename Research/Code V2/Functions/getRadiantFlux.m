function W = getRadiantFlux(Lm,lMin,lMax,lDelta,psd)

V = getEyeSens(lMin,lMax,lDelta,1978);
W=Lm/(683*sum(psd.*V)*lDelta);

end