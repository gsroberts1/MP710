cd('C:\Users\groberts\Desktop\301 - Slice 20 - L1');
load('stiff')
load('mask')
cd('Accel1p5')
load('stiffCS1p5');
stiff = stiff.*mask;
stiffCS1p5 = stiffCS1p5.*mask;
SSIM1p5 = ssim(stiffCS1p5,stiff)
save SSIM1p5

stiffNorm = stiff./(max(max(stiff)));
stiffNormCS1p5 = stiffCS1p5./(max(max(stiffCS1p5)));
PSNR1p5 = psnr(stiffNormCS1p5,stiffNorm)
save PSNR1p5

load('cols.mat');
load('rows.mat');
ROIstiff = stiff(rows,cols);
ROIstiffCS1p5 = stiffCS1p5(rows,cols);
ROI = [ROIstiff ROIstiffCS1p5];
n = size(ROIstiff,1)*size(ROIstiff,2);
sample1 = reshape(ROIstiff,1,n);
sample2 = reshape(ROIstiffCS1p5,1,n);
names = {'Fully Sampled Stiffness', 'Stiffness at R = 1.5', 'kPa'}; %%%
a = BlandAltman(sample1',sample2',names);