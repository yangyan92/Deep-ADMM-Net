function Gen_traindata( )
config;
ND = nnconfig.DataNmber;
%% Load samping pattern 
load('./mask/mask_20.mat')
save(strcat('./mask' , '.mat'), 'mask');
for i = 1:1:ND 
dir = './data/ChestTrain/im-';
load (strcat(dir , saveName(i,floor(log10(ND)))));
kspace_full = fft2(im_ori); 
y = (double(kspace_full)) .* (ifftshift(mask));
data.train = y;
data.label = im_ori;
save(strcat('./data/ChestTrain_sampling/', saveName(i, 2), '.mat'), 'data');
end