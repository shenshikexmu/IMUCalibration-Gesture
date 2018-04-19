


1.读入数据，calfata.mat

2.运行校正算法
   [Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm]=ImuCalibration_Gesture(cal_data)


3.校正部分：
   加速度、角速度，参照ICRA2014论文：<A Robust and Easy to implement method for imu
                                   calibration without External Equipments>
   磁力计部分，假设重力与磁向量的夹角不变。

4.参数部分：
  cal_acc=Ta*Ka*(raw_acc+Ba)
  cal_gyro=Tg*Kg*(raw_gyro+Bg)
  cal_mag=Tm2a*(raw_mag+Bm)
   
5.姿态部分：
   Mahony filter:
   conference <Nonlinear Complementery Filters on the Special Orthogonal Group>
   inspired by    http://blog.csdn.net/luoshi006/article/details/51513580
   EKF:
   derivation <A Double-Stage Kalman Filter for Orientation Tracking with 
               an Integrated Processor in 9-D IMU>
   high loe pass:
   high and low pass filter to Gyro atitude with Accelerate & Magnetic
   