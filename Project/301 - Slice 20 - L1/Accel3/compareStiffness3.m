cd('C:\Users\groberts\Desktop\301 - Slice 20 - L1');
load('stiff')
load('mask')
cd('Accel3') %%%
load('stiffCS3'); %%%
stiff = stiff.*mask;
stiffCS3 = stiffCS3.*mask; %%%
SSIM3 = ssim(stiffCS3,stiff) %%%
save SSIM3 %%%

stiffNorm = stiff./(max(max(stiff)));
stiffNormCS3 = stiffCS3./(max(max(stiffCS3))); %%%
PSNR3 = psnr(stiffNormCS3,stiffNorm) %%%
save PSNR3 %%%

load('cols.mat');
load('rows.mat');
ROIstiff = stiff(rows,cols);
ROIstiffCS3 = stiffCS3(rows,cols); %%%
ROI = [ROIstiff ROIstiffCS3]; %%%
n = size(ROIstiff,1)*size(ROIstiff,2);
sample1 = reshape(ROIstiff,1,n);
sample2 = reshape(ROIstiffCS3,1,n); %%%
names = {'Fully Sampled Stiffness', 'Stiffness at R = 3', 'kPa'}; %%%
a = BlandAltman(sample1',sample2',names);