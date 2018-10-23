% Plot varian power curve, given an array of reconstructed FIDs
%
% Samuel A. Hurley
% University of Wisconsin
% v1.1 26-Mar-2010
%
% Changelog:
%    v1.0 initial version (Dec-2009)
%    v1.1 include recon in plot_power.  Use subplots

function plot_power(fid)

if ~exist('fid', 'var')
  fid = 'fid';
end

pow = recon_1d(fid);

nsteps = size(pow, 3);
npts   = size(pow, 1);

for ii = 1:nsteps
  maxpow(ii) = max(squeeze(pow(:,1,ii)));
end

size(maxpow)

maxpow = circshift(maxpow', floor(nsteps/2));

pow = squeeze(pow);
pow = reshape(pow, [nsteps*npts 1]);
pow = circshift(pow, nsteps*npts/2);

figure;
subplot(2,1,1);
plot(pow);
title 'Prescan Power Steps'
subplot(2,1,2);
plot(maxpow);
title 'Prescan Power Curve'