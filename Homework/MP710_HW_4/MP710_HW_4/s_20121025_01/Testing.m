% % I. Problem One -- Reconstruct 2D IR Measurements and Generate T1 Maps %
% % The object of the first problem is to generate a T1 map of the phantom,
% % starting from the raw k-space data. % 1a.) Load Scan Information, Log,
% % and Header Into Matlab % Change directory into the series that contains
% % the inversion recovery T1 experiment. Display the 'text' file to verify
% % that we have used the correct pulse sequence. Display the scan log file
% % to verify that the experiment ran correctly without any warnings or
% % errors. % Then, load in the header and k-space data. Print the following
% % basic scan parameters to the command line: TR, TE, flip angle, and
% % inversion time (TI). Label these with the correct units.

clear all
close all

cd('../..')
% hw4dir = 'C:\Users\groberts\Documents\Fall 2018 Classwork\MP710\Homework\MP710_HW_4\MP710_HW_4\s_20121025_01';
hw4dir = 'C:\Users\robertsgr\Documents\MP710\Homework\MP710_HW_4\MP710_HW_4\s_20121025_01';
cd(hw4dir);                         % Start at home
addpath(genpath('.'))               % Add all subfolders to path
cd('./01.fid');                     % Go in the T1 mapping dataset

type('./text');                     % Pulse sequence info
type('./log');                      % Error log
info = load_procpar('./procpar');   % Header information

% Output TR, TE, flip angle, and TI:
disp(['The TR is: ' num2str(info.tr*1000) ' ms']);
disp(['The TE is: ' num2str(info.te*1000) ' ms']);
disp(['The flip angle is: ' num2str(info.flip1) ' degrees']);
ti = info.ti*1000;
disp(['The inversion times (TIs) are: ' num2str(floor(ti)) ' (ms)']);


% % 1b.) Load k-Space Data into Matlab and Reconstruct Images % Load the
% % raw k-space data file into Matlab using our custom command. Display the
% % matrix size on the command line. 

kspace = load_echoes('./fid');
data1 = fftshift(ifft2(fftshift(kspace))); 

mag1 = abs(data1);
phase1 = angle(data1);

figure;
subplot(1,2,1)
imshow(mag1(:,:,1),[])
title 'Magnitude';
subplot(1,2,2)
imshow(phase1(:,:,1),[])
title 'Phase';

% % 1c.) Investigate How the TI Parameter Affects Image Contrast & Signal %
% % In this section, we will investigate how the image intensity of each vial
% % changes with TI. From now on, we will be dealing only with magnitude
% % images.

maxPix = ceil(max(max(max(mag1))));
minPix = floor(min(min(min(mag1))));

figure;
for i=1:10
    subplot(2,5,i)
    imshow(mag1(:,:,i),[minPix maxPix])
end

% Vial #1 - Lower-Left   (0.40 mM Gd-DTPA)
figure;
imshow(mag1(:,:,1),[])
ROI1 = drawrectangle('Color','g','Label','ROI1');
title('Vial #1 - Lower-Left   (0.40 mM Gd-DTPA)');
i1 = floor(ROI1.Position);          % Position of the ROI [xmin, ymin, width, height]
                                    % xmin and ymin = upper left corner of ROI
rows1 = i1(2):(i1(2)+i1(4));
cols1 = i1(1):(i1(1)+i1(3));

means1 = zeros(1,10);
for i=1:10
    means1(i) = mean2(mag1(rows1,cols1,i));
end

sd1 = zeros(1,10);
for i=1:10
    sd1(i) = std2(mag1(rows1,cols1,i));
end

% Vial #2 - Upper-Middle (0.14 mM Gd-DTPA)
figure;
imshow(mag1(:,:,1),[])
ROI2 = drawrectangle('Color','g','Label','ROI2');
title('Vial #2 - Upper-Middle (0.14 mM Gd-DTPA)');
i2 = floor(ROI2.Position);          % Position of the ROI [xmin, ymin, width, height]
                                    % xmin and ymin = upper left corner of ROI
i2 = [i2(1) i2(2) i2(1)+i2(3) i2(2)+i2(4)];
rows2 = i2(2):i2(4);
cols2 = i2(1):i2(3);

means2 = zeros(1,10);
for i=1:10
    means2(i) = mean2(mag1(rows2,cols2,i));
end
sd2 = zeros(1,10);
for i=1:10
    sd2(i) = std2(mag1(rows2,cols2,i));
end

% Vial #3 - Lower-Right  (0.06 mM Gd-DTPA)
figure;
imshow(mag1(:,:,1),[])
ROI3 = drawrectangle('Color','g','Label','ROI3');
title('Vial #3 - Lower-Right  (0.06 mM Gd-DTPA)');
i3 = floor(ROI3.Position);          % Position of the ROI [xmin, ymin, width, height]
                                    % xmin and ymin = upper left corner of ROI
i3 = [i3(1) i3(2) i3(1)+i3(3) i3(2)+i3(4)];
rows3 = i3(2):i3(4);
cols3 = i3(1):i3(3);
means3 = zeros(1,10);
for i=1:10
    means3(i) = mean2(mag1(rows3,cols3,i));
end
sd3 = zeros(1,10);
for i=1:10
    sd3(i) = std2(mag1(rows3,cols3,i));
end

%Plot means of each ROI
figure;
subplot(3,1,1)
hold on
plot(ti,means1)
plot(ti,sd1)
title('Vial 1'); xlabel 'TI (ms)'; ylabel 'Signal [a.u.]'
hold off

subplot(3,1,2)
hold on
plot(ti,means2)
plot(ti,sd3)
title('Vial 2'); xlabel 'TI (ms)'; ylabel 'Signal [a.u.]'
hold off

subplot(3,1,3)
hold on
plot(ti,means3)
plot(ti,sd3)
title('Vial 3'); xlabel 'TI (ms)'; ylabel 'Signal [a.u.]'
hold off

% % 1d.) Fitting IR Signal to Compute T1 % Now we are going to fit the
% % inversion recovery data at each voxel to a mathematical model of the MRI
% % signal evolution. Instead of generating an MRI image of the phantom, with
% % different signal contrasts for each vial, we are now going to generate a
% % quantitative "map" showing estimates of the actual T1 times. The goal of
% % this question is to produce maps of proton density and T1. % Because we
% % have taken the magnitude operation on the data, we cannot simply fit it
% % to an exponential recovery curve, but instead must use nonlinear
% % least-squares fitting. % You will have to find a mathematical model for
% % the MRI signal as a function of inversion time (and any other relevant
% % parameters). These should be read directly from the image header. You
% % will also need to use a Matlab routine to fit the MRI data to this model,
% % and generate two images: PD and T1 HINT: If you are having trouble with
% % the fitting, try looking up Matlab documentation for "nonlinear
% % least-squares solver"
% Some code is provided to get started.

theta_init = [0 0];
pd = zeros(128,128);
t1 = zeros(128,128);
    
%% Loop over each voxel in the image.

for ii = 1:size(mag1, 1)
  for jj = 1:size(mag1, 2)
    
    % Keep track of progress
    progressbar(ii/(size(mag1,1)+1));
    
    % Grab the MRI data from each TI for this voxel
    vox_data = double(squeeze(abs(mag1(ii,jj,:))))';
    
    % Initial value TI value should be close to the proton density.
    theta_init(1) = mag1(ii,jj,1);
    % T1 should be approximately TI_null/0.69 when TR >> T1.
    theta_init(2) = min(vox_data)*1.443;

    % Signal model
    fun = @(theta) abs(theta(1).*(1-2*exp(-(ti./theta(2)))))-vox_data;
    
    % Finally, use a nonlinear least-squares solver to fit the data to the
    % model, using the initial guess as a starting point.
    lb = [0 10];
    ub = [30 6000];
    options = optimoptions('lsqnonlin','Display','off');
    theta = lsqnonlin(fun,theta_init,lb,ub,options);
    
    % Assign the results of the fitting to two output variables: pd and t1
   
    pd(ii,jj) = theta(1);
    t1(ii,jj) = theta(2);
  end
end

% Close the progress bar
progressbar(1);

% % 1e.) Analyze Results of T1 Mapping in a Variety of Ways Now that we have
% % fitted our MRI data to the inversion recovery model, we now want to
% % visualize and analyze our results in a number of ways. For starter's
% % let's look at images of PD and T1 maps. As before, display both images on
% % the same figure using subplot(). Put a color scale bar on each image and
% % label them with appropriate units. Mask out the background noise so that
% % it is easier to observe the actual phantom.

figure; imshow(pd,[]); colorbar; title('Proton Density Map');
figure; imshow(t1,[]); colorbar; title('T1 Map');

% Vial 1
pd_mean_vial1 = mean2(pd(rows1,cols1));
pd_sd_vial1 = std2(pd(rows1,cols1));
t1_mean_vial1 = mean2(t1(rows1,cols1));
t1_sd_vial1 = std2(t1(rows1,cols1));

% Vial 2
pd_mean_vial2 = mean2(pd(rows2,cols2));
pd_sd_vial2 = std2(pd(rows2,cols2));
t1_mean_vial2 = mean2(t1(rows2,cols2));
t1_sd_vial2 = std2(t1(rows2,cols2));

% Vial 3
pd_mean_vial3 = mean2(pd(rows3,cols3));
pd_sd_vial3 = std2(pd(rows3,cols3));
t1_mean_vial3 = mean2(t1(rows3,cols3));
t1_sd_vial3 = std2(t1(rows3,cols3));

Vial1 = [pd_mean_vial1; pd_sd_vial1; t1_mean_vial1; t1_sd_vial1];
Vial2 = [pd_mean_vial2; pd_sd_vial2; t1_mean_vial2; t1_sd_vial2];
Vial3 = [pd_mean_vial3; pd_sd_vial3; t1_mean_vial3; t1_sd_vial3];
Labels = {'PD Mean (a.u.)'; 'PD Standard Deviation (a.u.)'; 'T1 Mean (ms)'; 'T1 Standard Deviation (ms)'};

T = table(Vial1,Vial2,Vial3,'RowNames',Labels)

% Finally, display the fitted IR curve on top of the actual MRI data
% points (similar to what we did in 1c, but use the fitted results of PD
% and T1 to "fill in the curve" in-between actual data points).
% Label the axes with the proper names and units.

%% 
fit1 = abs(pd_mean_vial1.*(1-2*exp(-(ti./t1_mean_vial1))));
fit2 = abs(pd_mean_vial2.*(1-2*exp(-(ti./t1_mean_vial2))));
fit3 = abs(pd_mean_vial3.*(1-2*exp(-(ti./t1_mean_vial3))));

% Vial 1
figure;
hold on
plot(ti,fit1,'-.or');
scatter(ti,means1,'k','filled');
xlabel('TI (ms)'); ylabel('Signal (a.u.)'); title('Vial 1: Model vs Data'); legend('Fit', 'Empirical Data');
hold off

% Vial 2
figure;
hold on
plot(ti,fit2,'-.g');
scatter(ti,means2,'k','filled');
xlabel('TI (ms)'); ylabel('Signal (a.u.)'); title('Vial 2: Model vs Data'); legend('Fit', 'Empirical Data');
hold off

% Vial 3
figure;
hold on
plot(ti,fit3,'-.b');
scatter(ti,means3,'k','filled');
xlabel('TI (ms)'); ylabel('Signal (a.u.)'); title('Vial 3: Model vs Data'); legend('Fit', 'Empirical Data');
hold off

cd('../');
load_sdir;
%% 

% % II. Problem Two -- Reconstruct 2D Spine Echo Measurements and Generate T2
% % Maps The object of the second problem is to generate a T2 map of the
% % phantom, starting from the raw k-space data file. 2a.) Load Scan
% % Information, Log, and Header Into Matlab Change directory into the series
% % that contains the spin echo T2 mapping experiment. Display the 'text'
% % file to verify that we have used the correct pulse sequence. Display the
% % scan log file to verify that the experiment ran correctly without any
% % warnings or errors. Then, load in the header and k-space data. Print the
% % following basic scan parameters to the command line: TR, TEs, flip angle.
% % Label these with the correct units.

cd(hw4dir);                         % Start at home
cd('./02.fid');                     % Go in the T1 mapping dataset

type('./text');                     % Pulse sequence info
type('./log');                      % Error log
info2 = load_procpar('./procpar');  % Header information

% Output TR, TE, and flip angle
disp(['The TR is: ' num2str(info2.tr*1000) ' ms']);
te = info2.te*1000;
disp(['The TE times are: ' num2str(te) ' (ms)']);
disp(['The flip angle is: ' num2str(info2.flip1) ' degrees']);

%% 
% % 2b.) Load k-Space Data into Matlab and Reconstruct Images Load the raw
% % k-space data file into Matlab using our custom command. Display the
% % matrix size on the command line.

kspace2 = load_echoes('./fid');
disp(['The matrix size is: ' num2str(size(kspace2,1)) 'x' num2str(size(kspace2,2)) 'x' num2str(size(kspace2,3)) '.'])

% You will observe that the data size is 128x128x10. This corresponds to a
% readout length of 128 complex-valued points, 128 phase encode lines, and
% 10 different echo times (aka echoes, readouts, and blocks).
%
% As in (1b), perform a basic FFT reconstruction on the data to produce images.

data2 = fftshift(ifft2(fftshift(kspace2)));

mag2 = abs(data2);
phase2 = angle(data2);

% Next, let us take a look at our images to verify that our reconstruction looks
% good. Make a figure showing the magnitude and phase images of the first TE.
% Label the images them with appropriate titles and units.

figure;
subplot(1,2,1)
imshow(mag2(:,:,1),[])
title 'Magnitude (TE = 50 ms)';
subplot(1,2,2)
imshow(phase2(:,:,1),[])
title 'Phase (TE = 50 ms)';

% You should now have great looking magnitude and phase images.
%% 

% % 2c.) Investigate How the TE Parameter Affects Image Contrast & Signal In
% % this section, we will investigate how the image intensity of each vial
% % changes with echo time. From now on, we will be dealing only with
% % magnitude images.

% Start by plotting the image for each TE. Again, use subplot() so that all 10
% images are in the same figure, and label each of them with the appropriate TE.
% All images should be displayed on the same scale, and in gray-scale.

maxPix = ceil(max(max(max(mag2))));
minPix = floor(min(min(min(mag2))));

figure;
for i=1:10
    subplot(2,5,i)
    imshow(mag2(:,:,i),[minPix maxPix])
end

% Unlike the case for the IR experiment, increasing TE does not create "null
% points" or drastically change the contrast of signal between the vials.
% Instead, the overall signal intensity gradually declines as TE becomes longer.
%
% Let's make a plot showing the mean and standard deviation of MRI signal for
% each vial as a function of TE, so we can better visualize how the signal is
% changing. Try to put all three plots in the same figure if possible, and label
% the axes with the proper units. Use the ROIs defined from #1.

% Vial #1 - Lower-Left   (0.40 mM Gd-DTPA)
means1 = zeros(1,10);
for i=1:10
    means1(i) = mean2(mag2(rows1,cols1,i));
end
sd1 = zeros(1,10);
for i=1:10
    sd1(i) = std2(mag2(rows1,cols1,i));
end

% Vial #2 - Upper-Middle (0.14 mM Gd-DTPA)
means2 = zeros(1,10);
for i=1:10
    means2(i) = mean2(mag2(rows2,cols2,i));
end
sd2 = zeros(1,10);
for i=1:10
    sd2(i) = std2(mag2(rows2,cols2,i));
end

% Vial #3 - Lower-Right  (0.06 mM Gd-DTPA)
means3 = zeros(1,10);
for i=1:10
    means3(i) = mean2(mag2(rows3,cols3,i));
end
sd3 = zeros(1,10);
for i=1:10
    sd3(i) = std2(mag2(rows3,cols3,i));
end

%Plot means of each ROI
figure;
subplot(3,1,1)
hold on
plot(te,means1)
plot(te,sd1)
title('Vial 1'); xlabel 'TE (ms)'; ylabel 'Signal [a.u.]'
hold off

subplot(3,1,2)
hold on
plot(te,means2)
plot(te,sd3)
title('Vial 2'); xlabel 'TE (ms)'; ylabel 'Signal [a.u.]'
hold off

subplot(3,1,3)
hold on
plot(te,means3)
plot(te,sd3)
title('Vial 3'); xlabel 'TE (ms)'; ylabel 'Signal [a.u.]'
hold off

% As expected, the signal slowly decays away as TE increases. Although it may
% not be readily apparent, these data points follow the shape of a
% mono-exponential decay curve. We simply have not sampled a long enough TE to
% see the signal decay away to a value near zero.
%% 

% % 2d.) Fitting Spin-Echo Signal to Compute T2 Now we are going to fit the
% % spin echo data at each voxel to a mathematical model of the MRI signal
% % evolution. This time, the signal follows a more simple mathematical
% % function that can be cast into a linear form. Thus, we do not need to use
% % the complicated nonlinear least squares methods as done in Problem 1.
% % Instead, you should find a transform to make the signal linear, then use
% % standard least-squares regression to fit the model to the signal.
% As before, some code is provided to get started.

% Loop over each voxel in the image.

theta_init = [0 0];
pd = zeros(128,128);
t1 = zeros(128,128);

for ii = 1:size(mag2, 1)
  for jj = 1:size(mag2, 2)
    
    % Keep track of progress
    progressbar(ii/(size(mag2,1)+1));
    
    % Grab the MRI data from each TE for this voxel
    vox_data = double(squeeze(abs(mag2(ii,jj,:))))';
    
    % Initial value should be close to the proton density.
    theta_init(1) = mag2(ii,jj,1);
    % e^-x ~ 1 - x
    % Since M = M0*e^-(te/t2) we can use a first order approximation for
    % signals at M1 (M at te1) and M2 (M at te2) to approximate t2.
    % T2 ~ M1*TE2/(M1-M2)
    theta_init(2) = (mag2(ii,jj,1)*te(2))/(mag2(ii,jj,1)-mag2(ii,jj,2));

    % Signal model
    fun = @(theta) theta(1)*exp(-(te./theta(2)))-vox_data;
    
    % Finally, use a nonlinear least-squares solver to fit the data to the
    % model, using the initial guess as a starting point.
    lb = [0 0.001];
    ub = [30 3000];
    options = optimoptions('lsqnonlin','Display','off');
    theta = lsqnonlin(fun,theta_init,lb,ub,options);
    
    % Assign the results of the model to two output variables: pd and t2
    pd(ii,jj) = theta(1);
    t2(ii,jj) = theta(2);
  end
end

% Close the progress bar
progressbar(1);
%% 

% % 2e.) Analyze Results of T2 Mapping in a Variety of Ways Now that we have
% % fitted our MRI data to the spin echo model, we now want to visualize and
% % analyze our results in a number of ways. Display T2 and proton density
% % maps on the same figure using subplot(). Put a color scale bar on each
% % image and label them with appropriate units. Mask out the background
% % noise so it is easier to observe the actual phantom.

figure; subplot(1,2,1); imshow(pd,[]); colorbar; title('Proton Density Map');
subplot(1,2,2); imshow(t2,[]); colorbar; title('T2 Map');

% *The range of T2 should be between 0 and 600 milliseconds
% *Each of the three vials + background water should have a different T2 values.
%
% Compute the mean and standard deviation of both PD and T1 using the regions of
% interest from 1c. Display these values, making sure to label units!

% Vial 1
pd_mean_vial1 = mean2(pd(rows1,cols1));
pd_sd_vial1 = std2(pd(rows1,cols1));
t2_mean_vial1 = mean2(t2(rows1,cols1));
t2_sd_vial1 = std2(t2(rows1,cols1));

% Vial 2
pd_mean_vial2 = mean2(pd(rows2,cols2));
pd_sd_vial2 = std2(pd(rows2,cols2));
t2_mean_vial2 = mean2(t2(rows2,cols2));
t2_sd_vial2 = std2(t2(rows2,cols2));

% Vial 3
pd_mean_vial3 = mean2(pd(rows3,cols3));
pd_sd_vial3 = std2(pd(rows3,cols3));
t2_mean_vial3 = mean2(t2(rows3,cols3));
t2_sd_vial3 = std2(t2(rows3,cols3));

Vial1 = [pd_mean_vial1; pd_sd_vial1; t2_mean_vial1; t2_sd_vial1];
Vial2 = [pd_mean_vial2; pd_sd_vial2; t2_mean_vial2; t2_sd_vial2];
Vial3 = [pd_mean_vial3; pd_sd_vial3; t2_mean_vial3; t2_sd_vial3];
Labels = {'PD Mean (a.u.)'; 'PD Standard Deviation (a.u.)'; 'T2 Mean (ms)'; 'T2 Standard Deviation (ms)'};

T = table(Vial1,Vial2,Vial3,'RowNames',Labels)

%% 

% Finally, display the fitted spin echo signal on top of the actual MRI data
% points (similar to what we did in 1e, but for PD and T2).
% Label the ordinate and abscissa with the proper names and units.
% Extrapolate your fitted signal curve all the way out to 1000 ms.

fit1 = pd_mean_vial1*exp(-(te./t2_mean_vial1));
fit2 = pd_mean_vial2*exp(-(te./t2_mean_vial2));
fit3 = pd_mean_vial3*exp(-(te./t2_mean_vial3));

% Vial 1
figure;
hold on
plot(te,fit1,'-.or');
scatter(te,means1,'k','filled');
xlabel('T2 (ms)'); ylabel('Signal (a.u.)'); title('Vial 1: Model vs Data'); legend('Fit', 'Empirical Data');
hold off

% Vial 2
figure;
hold on
plot(te,fit2,'-.g');
scatter(te,means2,'k','filled');
xlabel('T2 (ms)'); ylabel('Signal (a.u.)'); title('Vial 2: Model vs Data'); legend('Fit', 'Empirical Data');
hold off

% Vial 3
figure;
hold on
plot(te,fit3,'-.b');
scatter(te,means3,'k','filled');
xlabel('T2 (ms)'); ylabel('Signal (a.u.)'); title('Vial 3: Model vs Data'); legend('Fit', 'Empirical Data');
hold off

% Congratulations! You have just performed your first T2 mapping experiment.
cd('..');

% % III - Appendix Series #1 and #2 are all that are needed to complete this
% % assignment. However, during the laboratory, we also acquired some
% % additional scans. These have been included in case you would like to play
% % around with the data. In particular, you should look at the k-space of
% % the EPI scan. It is very different from the scans we have just looked at.
% % EPI is difficult to reconstruct, however a custom command has been
% % provided to load the images that were reconstructed on the scanner:
% % epi_05 = load_fdf('05.img/slice001image001echo001.fdf',1); epi_06 =
% % load_fdf('06.img/slice001image001echo001.fdf',1); epi_07 =
% % load_fdf('07.img/slice001image001echo001.fdf',1);

figure;
subplot(1,3,1);
imagesc(epi_05);
axis image; axis off;
colormap gray;
title 'Gradient Echo EPI, 1-Shot';
subplot(1,3,2);
imagesc(epi_06);
axis image; axis off;
colormap gray;
title 'Gradient Echo EPI, 4-Shots';
subplot(1,3,3);
imagesc(epi_07);
axis image; axis off;
colormap gray;
title 'Spin Echo EPI, 4-Shots';

% % IV. -- Handing In The Assignment When you have competed this assignment
% % to your satisfaction, do the following: 1.) Save a copy of this M-file
% % with your name. 2.) Go to File->Publish. This will run the code over from
% % the very beginning, and will generate an html file of your code and png
% % images of your figures. Depending on your Matlab version, these may be
% % saved in a folder called 'html' 3.) Zip up these files, along with the
% % M-file, and e-mail them in to Professor Wieben