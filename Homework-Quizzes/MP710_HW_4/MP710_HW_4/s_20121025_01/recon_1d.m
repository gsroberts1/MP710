% FUNCTION [spect phase k1] = recon_1d(fid)
%
% Reconstruct a 1d free induction decay from a varian .fid file
%
% Samuel A. Hurley
% University of Wisconsin
% V1.2 26-Mar-2010
%
% Change History:
%    -1.1 updated documentation (Jul-2009)
%    -1.2 Automatically load fid if no arguments are supplied (Mar-2010)

function [spect phase k1] = recon_1d(fid)

% Specify default FID file
if ~exist('fid', 'var')
  fid = 'fid';
end

k1 = load_echoes(fid);
k1_fft = fftshift(fft(k1));

spect   = abs(k1_fft);
phase   = angle(k1_fft);

return;