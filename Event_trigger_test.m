clc
clear all
close all

aerfa=1;
x0=zeros(3,1);
x0=[96;4.8;0.2];
x1=zeros(3,1);
x1=[100;4.9;0.1];
x2=zeros(3,1);
x2=[90;4.7;0.1];
error1=x0-x1
error2=x0-x2-[10;0;0]
fai=[0,1.8,1.8];%增大fai，减小f
psai=[0,1.4,1.4];%触发条件中的误差权重向量,触发条件中的误差权重向量，增大pasi，增大f

%求解触发条件自适应参数
E=norm((psai.^0.5)*error1)*norm((psai.^0.5)*error1)
F=norm((fai.^0.5)*error2)*norm((fai.^0.5)*error2)
ri1=0.001;
ri2=0.0000005;
sgmai_min=1;
sgmai_max=2;%论文中是2
sgmai10=1;
sgmai20=1;

sgmai1=sgmai10/(1+ri1*sgmai10*E);
sgmai2=(sgmai20*E+ri2*sgmai_max)/(ri2+E);
sgma=sgmai1*aerfa+(1-aerfa)*sgmai2%求解触发条件自适应参数

%%触发条件
f=E-sgma*F
if f>0%%事件触发
    pp1=x1(1);
    vv1=x1(2);
    aa1=x1(3);
    x0=x1;
else
    pp1=x0(1);
    vv1=x0(2);
    aa1=x0(3);
end