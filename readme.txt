


1.读入数据，calfata.mat

2.运行校正算法
[Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm]=ImuCalibration_Gesture(cal_data)

3  校正部分：加速度、角速度，参照ICRA2014论文：A Robust and Easy to implement method for imu
   calibration without External Equipments
   磁力计部分，假设重力与磁向量的夹角不变。
   姿态部分：