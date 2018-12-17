cd('C:\Users\groberts\Desktop\301 - Slice 20 - L1');
load('stiff')
load('mask')
cd('Accel2')
load('stiffCS2');
stiff = stiff.*mask;
stiffCS2 = stiffCS2.*mask;
SSIM2 = ssim(stiffCS2,stiff)
save SSIM2

stiffNorm = stiff./(max(max(stiff)));
stiffNormCS2 = stiffCS2./(max(max(stiffCS2)));
PSNR2 = psnr(stiffNormCS2,stiffNorm)
save PSNR2

load('cols.mat');
load('rows.mat');
ROIstiff = stiff(rows,cols);
ROIstiffCS2 = stiffCS2(rows,cols);
ROI = [ROIstiff ROIstiffCS2];
n = size(ROIstiff,1)*size(ROIstiff,2);
sample1 = reshape(ROIstiff,1,n);
sample2 = reshape(ROIstiffCS2,1,n);
names = {'Fully Sampled Stiffness', 'Stiffness at R = 2', 'kPa'}; %%%
a = BlandAltman(sample1',sample2',names);