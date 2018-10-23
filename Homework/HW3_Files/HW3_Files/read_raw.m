function out = read_raw(in,Nx,Ny,im_no)
%function out = read_raw(in,Nx,Ny,im_no)
% Read raw k-space data from a General
%  Electric Signa raw data file. 
% This routine gets the data for one image
%  from a single or multi-slice data set
%
% Nx = readout length of the echo
% Ny = # of phase encodes
% im_no = # of image in sequence to be extracted
%
% written by Oliver Wieben, UW-Madison 21Mar2002
% 

% image number in the P-file
%  default = 1
if nargin <= 3,
  im_no = 1;
end

%% Important variables
 % byte size for each short
b_size = 2;
 % header size for LX file is 39984 bytes
header_size = 39984;

%% Calculate the offset in addition to the header
%  The scanner reserves space for calibration echoes
%   prior to every slice. These data are irrelevant 
%   since the calibration echoes are no longer acquired. 
offset = im_no*Nx*b_size*2 +(im_no-1)*Nx*Ny*b_size*2;

%% Go to the first data line of interest
%   and open big endian binary file
fid = fopen(in,'r','ieee-be');
if fid == -1
  disp(sprintf('File %s not found\n',fd));
  return;
end
fseek(fid, header_size + offset, -1);

%% Read the data in a temporary 1D array
%  The data are stored as short integers in 
%   inphase and quadrature pairs.
tmp1 = fread(fid,2*Nx*Ny, 'short');
fclose(fid);

% sort real and complex channels
tmp2 = tmp1(1:2:Nx*Ny*2) + i*tmp1(2:2:Nx*Ny*2);

%% Reshape into 2D matrix with Nx columns
%   and Ny rows
out = conj(reshape(tmp2,Nx,Ny)');





