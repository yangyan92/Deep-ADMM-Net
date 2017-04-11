# Deep-ADMM-Net

***********************************************************************************************************

This is a testing code for Basic-ADMM-Net in "Deep ADMM-Net for Compressive Sensing MRI" (NIPS 2016)
 
If you use this code, please cite our papers:
[1] Yan Yang, Jian Sun, Huibin Li, Zongben Xu. Deep ADMM-Net for Compressive Sensing MRI，NIPS(2016).

All rights are reserved by the authors.

Yan Yang -2017/04/05. More detail, feel free to contact: yangyan92@stu.xjtu.edu.cn


***********************************************************************************************************



Usage:


1. Load trained network with different stages in Basic_ADMM_Net_rec.m

   The models in './net/network_20' are trained from 100 real MR trainging images with 20% sampling rate. 
   The models in './net/network_30' are trained from 100 real MR trainging images with 30% sampling rate.
   
   If you apply these networks to  reconstruct of other MR images, it is best to re-train the models.

2. Load trained sampling pattern with different sampling ratios in Basic_ADMM_Net_rec.m

   The mask in './mask/mask_20' is a pseudo radial sampling patterns with 20% sampling rate.
   
3. Load test image  in Basic_ADMM_Net_rec.m

   The images in './data/Brain_data' are real-valued brain MR images.
   The images in './data/Chest_data' are 50 real-valued chest MR testing images in our paper.

4. Network setting is in  './config.m '.

5. To test our ADMM-Net, run
   ./Basic_ADMM_Net_rec.m





