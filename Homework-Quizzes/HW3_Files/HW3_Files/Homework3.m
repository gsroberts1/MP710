function MAG = get_mag(DATA)
%% GET_MAG Problem 1
MAG = sqrt((real(DATA)).^2 + (imag(DATA)).^2);
end
function PHASE = get_phase(DATA)
PHASE = atan(imag(DATA)./real(DATA));
end
%% 
% Part a)
%Got to HW3_Files directory
%Read in file and show image
file1 = [pwd '\test_pattern.png'];
data = importdata(file1);
figure; imshow(data); truesize
%Perform Fourier Trasnform (fast fourier)
DATA = fft2(data);
%Can we get back to the original image?
data_test = uint8(ifft2(DATA));
figure; imshow(data_test); truesize
test = isequal(data,data_test);
disp(test)
% If test = 1, then we got back the original image
realDATA = real(DATA);
imagDATA = imag(DATA);
MAG = sqrt(realDATA.^2 + imagDATA.^2);
imshow(log(MAG),[])
PHASE = atan(imagDATA./realDATA);
imshow(PHASE,[])  
%% 
% Part b)
[Nx,Ny] = size(DATA);
%How many lines do we need to exclude to reduce Ny by 50%?
index = (Ny*0.5)/2;
partial = [DATA(1:index,:); DATA((end-index+1):end,:)];
%Cut the data in half
% partial(1:index,:) = 0;
% partial((end-index+1):end,:)=0;
% %Looking at 1 column of data, samples should now be cut in half in the phase direction.
% disp(length(nonzeros(partial(:,1))))
% imshow(partial)