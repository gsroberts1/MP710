clear real imag
filter = window2(1600,1600,@gausswin);

% real = real(full);
% imag = imag(full);
% real = conv2(real,filter);
% imag = conv2(imag,filter);

NEW = padarray(FULL,[500 500],0);
FILT = filter.*NEW;
filt = ifft2(FILT);
imshow(filt)