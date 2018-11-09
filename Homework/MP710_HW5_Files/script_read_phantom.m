%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
% data path
% enter the proper directory here:
 path_root = 'E:\ZZZ  Med Phys 2014\05 Homeworks\HW3 - Phase\Posted\Data_Phantom\';

% the raw data files stored in Siemens Sonata format 2004
im_path(1,:) = [path_root '01\meas.out'];
im_path(2,:) = [path_root '02\meas.out'];
im_path(3,:) = [path_root '03\meas.out'];

% debug level: 0 <= debug <= 5 
debug = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF INPUTS

%%% read all measured data into matrix raw_data
[raw_vol,im_vol] = read_from_raw_n4(im_path,debug);




