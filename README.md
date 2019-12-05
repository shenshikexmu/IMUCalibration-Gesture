# IMUCalibration-Gesture
calibration for Imu and show gesture
#### 0.相关博客:

https://blog.csdn.net/shenshikexmu/article/details/80013444

#### 1.读入数据

 load('calfata.mat')

#### 2.运行校正算法

   [Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm,mag_strength]=ImuCalibration_Gesture(cal_data)


#### 3.校正部分

##### 加速度、角速度
   conference  **A Robust and Easy to implement method for imu calibration without External Equipments**
##### 磁力计
   算法mag2acc_matrix假设重力与磁向量的夹角不变，算法Cal_mag4acc_frame利用不同姿态下传感器感受的磁通向量的变化与姿态变化的相关性，计算参数。

#### 4.参数部分

#####  cal_acc=Ta*Ka*(raw_acc+Ba)
#####  cal_gyro=Tg*Kg*(raw_gyro+Bg)
#####  cal_mag=Tm2a*(raw_mag+Bm)
   
#### 5.姿态部分

 #####  Mahony filter
   conference **Nonlinear Complementery Filters on the Special Orthogonal Group**
   
   inspired by    http://blog.csdn.net/luoshi006/article/details/51513580
 #####  EKF
   derivation **A Double-Stage Kalman Filter for Orientation Tracking with an Integrated Processor in 9-D IMU**
 #####  high low pass
   high and low pass filter to Gyro atitude with Accelerate & Magnetic
 ##### 滤波后的四元数
 
![滤波结果](https://img-blog.csdn.net/20180606120722833?watermark/2/text/aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3NoZW5zaGlrZXhtdQ==/font/5a6L5L2T/fontsize/400/fill/I0JBQkFCMA==/dissolve/10)
##### 算法视频
[![视频]( https://vthumb.youku.com/054208085DE892CC000001630908ED3D)](https://v.youku.com/v_show/id_XNDQ1ODc2MzAyNA==.html?sharefrom=iphone&sharekey=375f30a8f7e787f23d9927c39e6d9be38)

