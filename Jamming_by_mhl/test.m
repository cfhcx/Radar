clear;
close all;
clc;

PRI=100e-6;
% PRI=50e-6;
DutyRatio=0.3;
Tp=PRI*DutyRatio;
JammingPower=2;%��������������
CarrierPower=30;%��������������
NoisePower=28;%��������������������������׼�� ��������Ϊ����
sigma=2;%���������� �����ı�׼��
R=8e3;                            
c=3e8;

B=15e6;%�״��źŴ���
fs=2*B;
Bj=10e6;%�����źŴ���
fj=15e6;%��������
f=40e6;%�״��ź�
% Kfm=1.66*sigma*Bj;%��Ƶб��
Kfm=7*1e6;
Kpm=1e-1;%����б��
D=Kpm*sigma;

PulseNum=5; %SMSP��
RectPulseNum=9;
TimeSlotNum=3;
% Tp=PRP*DutyRatio;
% N_Tp=fix(Tp*fs);
% t_Tp=(0:N_Tp-1)/fs;
T1=1*PRI;T2=25*PRI;T3=1*PRI;T4=1*PRI;
% T1=0.001;T2=0.1;T3=0.001;T4=0.001;
v_f=2e7; %�ٶȲ��������ٶ�m/s
v_t=1e6; %���벨�������ٶ�m/s
v=1e3; %m/s
Bd=5e5;%�����յ��ƴ��� ������ɱȽϴ��ֵ ���ӻ���ͼ��Խ�����Եķ�Ӧ��������Ƶ��
fd0=2e5;%�ٶ�����Ƶ�ʵ���ʼֵ
% Bd=5;%�����յ��ƴ���  �����Ķ�����Ƶ�ʺ�С
% fd0=2;%�ٶ�����Ƶ�ʵ���ʼֵ
delta_R=1e2;
n=7;%��Ŀ�����
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
a0=0.5;%��������ռ�ձ�
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ʦ�ֵĴ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wr=4e3;
% T=Wr/(c/2);
% nrn=fix((T+Tp)*fs);
% prf=1/PRI;
% lambda=c/f;
% resolution_range=c/2/B;
% Da=2*resolution_range;                      %%%���߿׾�����
% Cta_a=lambda/Da*180/pi;                     %%%������Ƚ�
% Ts=Cta_a*T/lambda;                          %%%�ϳɿ׾�ʱ��
% Bd=2*v*Cta_a/lambda;                        %%%�����յ��ƴ���
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
% title('ʱƵͼ','fontsize',12,'fontweight','bold');
% xlabel('ʱ��(us)','fontsize',13,'fontweight','bold');
% ylabel('Ƶ��(MHz)','fontsize',13,'fontweight','bold');

t=(0:length(J)-1)/fs;
figure;plot(t*1e6,real(J));
title('�����ź�ʱ��ͼ','fontsize',12,'fontweight','bold');
xlabel('ʱ��(us)','fontsize',12,'fontweight','bold');
ylabel('����','fontsize',12,'fontweight','bold');

% F=fftshift(fft(J));
F=fft(J);
freq=(0:length(F)-1)*fs/length(F);
figure;plot(freq*1e-6,abs(F));
title('�����ź�Ƶ��ͼ','fontsize',12,'fontweight','bold');
xlabel('Ƶ��(MHz)','fontsize',12,'fontweight','bold');
ylabel('����','fontsize',12,'fontweight','bold');


