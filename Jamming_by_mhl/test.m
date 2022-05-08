clear;
close all;
clc;

PRI=100e-6;
% PRI=50e-6;
DutyRatio=0.3;
Tp=PRI*DutyRatio;
JammingPower=2;%仿真中做幅度用
CarrierPower=30;%仿真中做幅度用
NoisePower=28;%噪声调幅仿真中做幅度用作标准差 噪声功率为方差
sigma=2;%噪声调相用 噪声的标准差
R=8e3;                            
c=3e8;

B=15e6;%雷达信号带宽
fs=2*B;
Bj=10e6;%干扰信号带宽
fj=15e6;%调制噪声
f=40e6;%雷达信号
% Kfm=1.66*sigma*Bj;%调频斜率
Kfm=7*1e6;
Kpm=1e-1;%调相斜率
D=Kpm*sigma;

PulseNum=5; %SMSP用
RectPulseNum=9;
TimeSlotNum=3;
% Tp=PRP*DutyRatio;
% N_Tp=fix(Tp*fs);
% t_Tp=(0:N_Tp-1)/fs;
T1=1*PRI;T2=25*PRI;T3=1*PRI;T4=1*PRI;
% T1=0.001;T2=0.1;T3=0.001;T4=0.001;
v_f=2e7; %速度波门拖引速度m/s
v_t=1e6; %距离波门拖引速度m/s
v=1e3; %m/s
Bd=5e5;%多普勒调制带宽 可以设成比较大的值 可视化的图像越能明显的反应出多普勒频移
fd0=2e5;%假多普勒频率的起始值
% Bd=5;%多普勒调制带宽  正常的多普勒频率很小
% fd0=2;%假多普勒频率的起始值
delta_R=1e2;
n=7;%假目标个数
% AM_jamming=AM_jamming(PRI,DutyRatio,CarrierPower,NoisePower,Bj,fs,fj,R);
% J=AM_jamming;
% FM_jamming=FM_jamming(PRI,DutyRatio,JammingPower,sigma,Kfm,Bj,fs,fj,R);
% J=FM_jamming;
% PM_jamming=PM_jamming(PRI,DutyRatio,JammingPower,sigma,Kpm,Bj,fs,fj,R);
% J=PM_jamming;
% SMSP_jamming=SMSP_jamming(PRI,DutyRatio,JammingPower,B,fs,f,PulseNum,R);
% J=SMSP_jamming;
% CandI_jamming=CandI_jamming(PRI,DutyRatio,JammingPower,B,fs,fj,RectPulseNum,TimeSlotNum,R);
% J=CandI_jamming;
% VGP_jamming=VGP_jamming(PRI,DutyRatio,JammingPower,B,fs,fj,T1,T2,T3,T4,v_f,v,R);
% J=VGP_jamming;
% RGP_jamming=RGP_jamming(PRI,DutyRatio,JammingPower,B,fs,fj,T1,T2,T3,T4,v_t,v,R);
% J=RGP_jamming;
delta_t=Tp/20;
a0=0.5;%矩形脉冲占空比
% TDSCR_jamming=TDSCT_jamming(PRI,DutyRatio,JammingPower,B,fs,f,delta_t,R);
% J=TDSCR_jamming;
% TDSCR1_jamming=TDSCT1_jamming(PRI,DutyRatio,JammingPower,B,fs,f,delta_t,a0,R);
% J=TDSCR1_jamming;
% NC_jamming=NC_jamming(PRI,DutyRatio,JammingPower,B,Bj,fs,f,R);
% J=NC_jamming;
% MFT_jamming=MFT_jamming(PRI,DutyRatio,JammingPower,B,Bd,fd0,fs,f,delta_R,R,n);
% J=MFT_jamming;
JSR=10; radarPower=1; 
MFT_V_jamming=MFT_V_jamming(PRI,DutyRatio,JSR,radarPower,B,Bd,fd0,fs,f,delta_R,R,n);
J=MFT_V_jamming;
% PE_jamming=PE_jamming(PRI,DutyRatio,JammingPower,B,fs,f);
% J=PE_jamming;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%师兄的代码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wr=4e3;
% T=Wr/(c/2);
% nrn=fix((T+Tp)*fs);
% prf=1/PRI;
% lambda=c/f;
% resolution_range=c/2/B;
% Da=2*resolution_range;                      %%%天线孔径长度
% Cta_a=lambda/Da*180/pi;                     %%%波束宽度角
% Ts=Cta_a*T/lambda;                          %%%合成孔径时间
% Bd=2*v*Cta_a/lambda;                        %%%多普勒调制带宽
% % Ka=2*v^2/lambda/R;  Ka=Bd/Ts;                        
% nan_num=prf*lambda*R/(2*v*resolution_range);
% data_len=fix(PRI*fs);
% SF_jamming=SF_jamming(Tp,B,fs,nrn,nan_num,data_len);
% J=SF_jamming;
% SWSF_jamming=SWSF_jamming(Tp,B,fs,nrn,nan_num,data_len);
% J=SWSF_jamming;
% PSF_jamming=PSF_jamming(Tp,B,fs,nrn,nan_num,data_len);
% J=PSF_jamming;
% ISR_jamming=ISR_jamming(Tp,B,fs,nrn,nan_num,data_len);
% J=ISR_jamming;
% RT_jamming=RT_jamming(Tp,B,fs,nrn,nan_num,data_len);
% J=RT_jamming;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [tfr,t,f]=tfrstft(J',1:length(J),1024); 
% % figure;imagesc(t/fs*1e6,-f(1:length(f)/2+1)*2*fs*1e-6,abs(tfr));
% figure;imagesc(t/fs*1e6,-f/2*fs,abs(tfr));
% % ylim([0 1024*fs*1e-3]); 
% axis xy
% % set(gca,'YTick',0:fs*1e-3:1024*fs*1e-3);
% title('时频图','fontsize',12,'fontweight','bold');
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('频率(MHz)','fontsize',13,'fontweight','bold');

t=(0:length(J)-1)/fs;
figure;plot(t*1e6,real(J));
title('干扰信号时域图','fontsize',12,'fontweight','bold');
xlabel('时间(us)','fontsize',12,'fontweight','bold');
ylabel('幅度','fontsize',12,'fontweight','bold');

% F=fftshift(fft(J));
F=fft(J);
freq=(0:length(F)-1)*fs/length(F);
figure;plot(freq*1e-6,abs(F));
title('干扰信号频谱图','fontsize',12,'fontweight','bold');
xlabel('频率(MHz)','fontsize',12,'fontweight','bold');
ylabel('幅度','fontsize',12,'fontweight','bold');


