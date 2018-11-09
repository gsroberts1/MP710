function [raw,im] = read_from_raw_n4(data_path,debug)
%function [raw,im] = read_from_raw_n4_local(in,N_read,debug)
% Read raw data from a Numaris4 scan (VA21 or 23) into a complex matrix.  
%  The function takes the binary file IN (a Numaris4 file 'meas.out'), 
%  strips off all the header information, and stores the data into a 
%  a complex 2D matrix OUT.
% In the case of multi-receiver acquisitions, the data are stored 'interleaved'.
%  For example with four receiver channels the first echo is stored four
%  times before the second echo is stored:
%    echo 1, receiver 1
%    echo 1, receiver 2
%    echo 1, receiver 3
%    echo 1, receiver 4
% No test on 3D sequences yet!
%
% DEBUG defines the level of debug information during execution.
%
% Example: out = read_raw_n4('meas.out',5);
%
% written by Oliver Wieben, Uniklinik Freiburg 12Sep2002
%      10Nov03 OW automatic detection of # of samples per readout
% 

no_im=size(data_path,1);

%%% read all measured data into matrix raw_data_full
for i= 1:no_im
    raw_full(:,:,i)=read_raw_n4_local(data_path(i,:),0);    
end

%%% downsample the array 
Nx=size(raw_full,2);
raw(:,:,:)=raw_full(:,1:2:Nx,:);

%% calculate the 2D inverse FT for echo inspection
for i= 1:no_im
    temp = squeeze(raw(:,:,i));
    im(:,:,i) = fftshift(fft2(fftshift(temp)));
end

if debug > 0
  fprintf('\nRead image volume: %d x %d x %d\n',size(raw,1),size(raw,2),size(raw,3))    
end


function out = read_raw_n4_local(in,debug)
%function out = read_raw_n4_local(in,N_read,debug)
% Read raw data from a Numaris4 scan (VA21 or 23) into a complex matrix.  
%  The function takes the binary file IN (a Numaris4 file 'meas.out'), 
%  strips off all the header information, and stores the data into a 
%  a complex 2D matrix OUT.
% In the case of multi-receiver acquisitions, the data are stored 'interleaved'.
%  For example with four receiver channels the first echo is stored four
%  times before the second echo is stored:
%    echo 1, receiver 1
%    echo 1, receiver 2
%    echo 1, receiver 3
%    echo 1, receiver 4
% No test on 3D sequences yet!
%
% DEBUG defines the level of debug information during execution.
%
% Example: out = read_raw_n4('meas.out',5);
%
% written by Oliver Wieben, Uniklinik Freiburg 12Sep2002
%      10Nov03 OW automatic detection of # of samples per readout
% 

if debug == -1
  in = 'd:\user\wieben\Work_in_Progress\Radiale_Bildgebung\Daten\Daten_12Sep02';
end
  
%% Define important variables
 % byte size per float = 4
b_float = 4;
%%% There are three types of additional headers information in the raw data files
 % 1. once at the beginning of the raw data file: header with 32 bytes
 %     contents: ???
h_raw_start = 32;
 % 2. in front of every echo: 128 bytes 
 %    contents: you can see the contens with mdb.exe, includes free parameters
h_read = 128;
 % 3. once at the end of the file: 384 bytes
 %    contens: ???
h_raw_end = 384;

%%%%%%%%%%%%%%%%% Step 1: Read all raw data from file and strip headers from start and end of file
fid = fopen(in,'r');
if fid == -1
  disp(sprintf('File %s not found\n',in));
  return;
end

%%%% Task 1: find the number of points per readout
 % go to first data line
offset_line = 0;
 % goto byte offset for entry 'samples in scan' 
offset_h = 28;
fseek(fid, h_raw_start + offset_line + offset_h, -1);
N_read = fread(fid, 1,'uint16');
fclose(fid);

%%% Not the greatest programming here: reopen the file again
fid = fopen(in,'r');
% skip over the first short header
fseek(fid, h_raw_start, -1);
 % read all data from file
[tmp,count] = fread(fid, 'float');
fclose(fid);
 % remove the bytes (header) from the end of the file
tmp = tmp(1:count-h_raw_end/b_float);

%%%%%%%%%%%%%%%%% Step 2: Remove headers from each echo 
 % # of floats read in
ld = length(tmp);
 % get the # of echoes: 
 % consider the bytes used for the headers of each echo in 'units' float
h_read_f = h_read/b_float;
N_e = floor(ld/(N_read*2+h_read_f));

 %% format the 1D float array into 2D complex matrix
tmp = tmp(1:2:N_e*(N_read*2+h_read_f)) + j*tmp(2:2:N_e*(N_read*2+h_read_f));
out = conj(reshape(tmp,N_read+h_read_f/2,N_e)');
 
 % finally: removal of the headers in front of every echo
out = out(:,h_read_f/2+1:N_read+h_read_f/2);

if debug > 0
  fprintf('\n Data size: N_read = %d, N_phase = %d',N_read,size(out,1));
end
  
%%% Example for output
if debug >  3
 imagesc(log(abs(out)))
 colorbar
 colormap('gray')
end
