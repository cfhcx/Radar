clear;
close all;
clc;

CarrierPower=30;NoisePower=10;
PRI=100e-6;DutyRatio=0.2;
Bj=10e6;B=10e6;fj=15e6;f=15e6;fs=30e6; R=8e3;

AM_jamming=AM_jamming(PRI,DutyRatio,CarrierPower,NoisePower,Bj,fs,fj,R);
J=AM_jamming;

t=(0:length(J)-1)/fs;
subplot(1,2,1);plot(t*1e6,real(J));
title('干扰信号时域图','fontsize',12,'fontweight','bold');
xlabel('时间(us)','fontsize',12,'fontweight','bold');
ylabel('幅度','fontsize',12,'fontweight','bold');

F=fft(J);
freq=(0:length(F)-1)*fs/length(F);
subplot(1,2,2);plot(freq*1e-6,abs(F));
title('干扰信号频谱图','fontsize',12,'fontweight','bold');
xlabel('频率(MHz)','fontsize',12,'fontweight','bold');
ylabel('幅度','fontsize',12,'fontweight','bold');