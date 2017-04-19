# Deep-ADMM-Net

***********************************************************************************************************

This is a testing and training code for Deep ADMM-Net in "Deep ADMM-Net for Compressive Sensing MRI" (NIPS 2016)<br>
 
If you use this code, please cite our paper:<br>
[1] Yan Yang, Jian Sun, Huibin Li, Zongben Xu. Deep ADMM-Net for Compressive Sensing MRI锛孨IPS(2016).<br>
[http://gr.xjtu.edu.cn/web/jiansun/publications](http://gr.xjtu.edu.cn/web/jiansun/publications) <br>
All rights are reserved by the authors.<br>

Yan Yang -2017/04/05. For more detail, feel free to contact: yangyan92@stu.xjtu.edu.cn<br>


***********************************************************************************************************



Usage:<br>

1.For testing the trained network <br>

1). Load trained network with different stages in 'main_ADMM_Net_test.m'

   The models in './net/network_20' are trained from 100 real MR trainging images with 20% sampling rate. 

<br>
   The models in './net/network_30' are trained from 100 real MR trainging images with 30% sampling 

rate.<br>
   
   If you apply ADMM-Net to  reconstruct  other MR images, it is best to re-train the models.<br>

2). Load  sampling pattern with different sampling ratios in  'main_ADMM_Net_test.m'<br>

   The mask in './mask/mask_20' is a pseudo radial sampling pattern with 20% sampling rate.<br>
   
3). Load test image  in  'main_ADMM_Net_test.m'.

   The images in './data/Brain_data' are real-valued brain MR images.
   The images in './data/Chest_data' are 50 real-valued chest MR testing images in our paper.

4). Network setting is in  'config.m '.

5). To test our ADMM-Net, run 'main_ADMM_Net_test.m'


<br>
<br>
2.For training the networks<br>
1). The training chest dataset is in './data/ChestTrain_20'.<br>
    Run "Gen_traindata.m" to generate training data, <br>
    and load  corresponding sampling pattern in this operation. <br>


2). Modify the network setting and trainging setting in  'config.m '.

3). To train ADMM-Net by L-BFGS algorithm, run 'main_netTrain.m' .

4). After training, the trained network and the training error are saved in './Train_output'.<br>



***********************************************************************************************************
 
The testing result of the demo images.

1) Brain_data1.（20% sampling rate）

|--------------|  re_LOss  |  re_PSnr  |  <br>
|  net-stage7- |  0.0578   |  35.60    |  <br>
|  net-stage14 |  0.0562   |  35.83    |  <br>
|  net-stage15 |  0.0561   |  35.85    |  <br>


2) Brain_data2.(20% sampling rate）

|--------------|  re_LOss  |  re_PSnr  |  <br>
|  net-stage7- |  0.0957   |  30.40    |  <br>
|  net-stage14 |  0.0929   |  30.65    |  <br>
|  net-stage15 |  0.0927   |  30.67    |  <br>










