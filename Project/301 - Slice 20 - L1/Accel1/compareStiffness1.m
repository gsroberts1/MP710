cd('C:\Users\groberts\Desktop\301 - Slice 20 - L1');
load('stiff')
load('mask')
cd('Accel1')
load('stiffCS1');
stiff = stiff.*mask;
stiffCS1 = stiffCS1.*mask;
SSIM = ssim(stiffCS1,stiff)
save SSIM1

stiffNorm = stiff./(max(max(stiff)));
stiffNormCS1p5 = stiffCS1./(max(max(stiffCS1)));
PSNR1 = psnr(stiffNormCS1p5,stiffNorm)
save PSNR1

load('cols.mat');
load('rows.mat');
ROIstiff = stiff(rows,cols);
ROIstiffCS1 = stiffCS1(rows,cols);
ROI = [ROIstiff ROIstiffCS1];
n = size(ROIstiff,1)*size(ROIstiff,2);
sample1 = reshape(ROIstiff,1,n);
sample2 = reshape(ROIstiffCS1,1,n);
names = {'Fully Sampled Stiffness', 'Stiffness at R = 1', 'kPa'}; %%%
a = BlandAltman(sample1',sample2',names);