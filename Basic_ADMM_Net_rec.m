function  Basic_ADMM_Net_rec(  )

% This is a cpu test code demo for ADMM_Net_v1 reconstruction.
% Output: the average NMSE and PSNR over the test images.

% clc;
% clear all;

%% Load trained network
load('./net/network_20/net-stage15.mat')
%% Load data 
load('./data/Brain_data/Brain_data2.mat')
%load('.data/Chest_data/chest_data1.mat')
load('./mask/mask_20.mat')
%% Undersampling in the k-space
kspace_full = fft2(im_ori); 
y = (double(kspace_full)) .* (ifftshift(mask));
data.train = y;
data.label = im_ori;

%% reconstrction by ADMM-Net
%tic
[re_LOss, rec_image] = loss_with_gradient_single_before(data, net);
%Time_Net_rec = toc
re_PSnr = psnr(abs(rec_image) , abs(data.label))
re_LOss
Zero_filling_rec = ifft2(y);
figure;
subplot(1,2,1); imshow(abs(Zero_filling_rec)); xlabel('Zero-filling reconstructon result');
subplot(1,2,2); imshow(abs(rec_image)); xlabel('ADMM-Net reconstruction result');
imwrite(abs(rec_image),'rec_image.png')
end

