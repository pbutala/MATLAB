function Lm = getLuminousFlux(W,lMin,lMax,lDelta,psd)

V = getEyeSens(lMin,lMax,lDelta,1978);
Lm=W*(683*sum(psd.*V)*lDelta);

end