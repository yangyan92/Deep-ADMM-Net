ï»? Deep-ADMM-Net

***********************************************************************************************************

This is a testing and training code for Basic-ADMM-Net in "Deep ADMM-Net for Compressive Sensing MRI" (NIPS 2016)
 
If you use this code, please cite our papers:
[1] Yan Yang, Jian Sun, Huibin Li, Zongben Xu. Deep ADMM-Net for Compressive Sensing MRIï¼ŒNIPS(2016).

All rights are reserved by the authors.

Yan Yang -2017/04/05. More detail, feel free to contact: yangyan92@stu.xjtu.edu.cn


***********************************************************************************************************



Usage:

1.For testing the trained network 

1). Load trained network with different stages in main_ADMM_Net_test.m

   The models in './net/network_20' are trained from 100 real MR trainging images with 20% sampling rate. 
   The models in './net/network_30' are trained from 100 real MR trainging images with 30% sampling rate.
   
   If you apply these networks to  reconstruct of other MR images, it is best to re-train the models.

2). Load trained sampling pattern with different sampling ratios in main_ADMM_Net_test.m

   The mask in './mask/mask_20' is a pseudo radial sampling patterns with 20% sampling rate.
   
3). Load test image  in main_ADMM_Net_test.m

   The images in './data/Brain_data' are real-valued brain MR images.
   The images in './data/Chest_data' are 50 real-valued chest MR testing images in our paper.

4). Network setting is in  'config.m '.

5). To test our ADMM-Net, run 'main_ADMM_Net_test.m'



2.For training the networks

1). The training chest dataset is in './data/ChestTrain_20'.
    Run 
"Gen_traindata.m" to generate training data,
    and load  corresponding sampling pattern in these operation. 

2). Modify the network setting and trainging setting in  'config.m '.

3). To train ADMM-Net by L-BFGS algorithm, run
 'main_netTrain.m' .

4). After training, the trained network and the training error are saved in './Train_output'.



***********************************************************************************************************
 
The testing result of the demo images.

1) Brain_data1.ï¼?0% sampling rateï¼?

| ------------ |  re_LOss  |  re_PSnr  |
|  net-stage7  |  0.0578   |  35.60    |
|  net-stage14 |  0.0562   |  35.83    |
|  net-stage15 |  0.0561   |  35.85    |


2) Brain_data2.(20% sampling rateï¼?

| ------------ |  re_LOss  |  re_PSnr  |
|  net-stage7  |  0.0957   |  30.40    |
|  net-stage14 |  0.0929   |  30.65    |
|  net-stage15 |  0.0927   |  30.67    |





