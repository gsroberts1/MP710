filter = window2(600,600,@hanning);

FILT = filter.*FULL;
filt = ifft2(FILT);
figure; imshow(abs(filt),[]);
figure; imshow(abs(full),[]);

