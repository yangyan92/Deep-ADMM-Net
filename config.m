%% network setting
global nnconfig;
nnconfig.FilterNumber =8;
nnconfig.FilterSize = 3 ;
nnconfig.Stage = 5;
nnconfig.Padding = 1;
nnconfig.LinearLabel = double(-1:0.02:1);
%% training and testing setting
nnconfig.EnableGPU = 0;    %first run on GPU always takes longer than subsequent runs
nnconfig.WeightDecay = 0;




%% data
nnconfig.TrainNumber = 1;
nnconfig.ValDataNumber = 50;
nnconfig.TesDataNumber = 1;
nnconfig.ImageSize = [256,256];
nnconfig.DataNmber = 100;
 




