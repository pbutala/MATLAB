function L = getLorentzian(MN,FWHM,X)
HWHM = FWHM/2;
L = HWHM*(1./((X-MN).^2 + HWHM^2))/pi;
end