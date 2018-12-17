SSIM1 = zeros(2,8);
SSIM1p5 = zeros(2,8);
SSIM2 = zeros(2,8);
SSIM3 = zeros(2,8);
for i = 1:2
    for j = 1:8
        if i == 1
            mag = mag_dat(:,:,j);
            SSIM1(i,j) = ssim(magCS1(:,:,j),mag);
            SSIM1p5(i,j) = ssim(magCS1p5(:,:,j),mag);
            SSIM2(i,j) = ssim(magCS2(:,:,j),mag);
            SSIM3(i,j) = ssim(magCS3(:,:,j),mag);
        else 
            phase = phase_dat(:,:,j);
            SSIM1(i,j) = ssim(phaseCS1(:,:,j),phase);
            SSIM1p5(i,j) = ssim(phaseCS1p5(:,:,j),phase);
            SSIM2(i,j) = ssim(phaseCS2(:,:,j),phase);
            SSIM3(i,j) = ssim(phaseCS3(:,:,j),phase);
        end 
    end 
end 

PSNR1 = zeros(2,8);
PSNR1p5 = zeros(2,8);
PSNR2 = zeros(2,8);
PSNR3 = zeros(2,8);
for i = 1:2
    for j = 1:8
        if i == 1
            mag = mag_dat(:,:,j);
            PSNR1(i,j) = psnr(magCS1(:,:,j),mag);
            PSNR1p5(i,j) = psnr(magCS1p5(:,:,j),mag);
            PSNR2(i,j) = psnr(magCS2(:,:,j),mag);
            PSNR3(i,j) = psnr(magCS3(:,:,j),mag);
        else 
            phase = phase_dat(:,:,j);
            PSNR1(i,j) = psnr(phaseCS1(:,:,j),phase);
            PSNR1p5(i,j) = psnr(phaseCS1p5(:,:,j),phase);
            PSNR2(i,j) = psnr(phaseCS2(:,:,j),phase);
            PSNR3(i,j) = psnr(phaseCS3(:,:,j),phase);
        end 
    end 
end 

mSSIM1 = mean(SSIM1,2)
mSSIM1p5 = mean(SSIM1p5,2)
mSSIM2 = mean(SSIM2,2)
mSSIM3 = mean(SSIM3,2)
mPSNR1 = mean(PSNR1,2)
mPSNR1p5 = mean(PSNR1p5,2)
mPSNR2 = mean(PSNR2,2)
mPSNR3 = mean(PSNR3,2)

ssimEPImag = [SSIM1(1) SSIM1p5(1) SSIM2(1) SSIM3(1)];
ssimEPIphase = [SSIM1(2) SSIM1p5(2) SSIM2(2) SSIM3(2)];
psnrEPImag = [PSNR1(1) PSNR1p5(1) PSNR2(1) PSNR3(1)];
psnrEPIphase = [PSNR1(2) PSNR1p5(2) PSNR2(2) PSNR3(2)];

