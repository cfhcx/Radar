function FM_jamming=FM_jamming(PRI,DutyRatio,JammingPower,sigma,Kfm,Bj,fs,fj,R)
%PRI脉冲重复周期；DutyRatio占空比；JammingPower干扰功率；Bj干扰带宽；
%fs采样频率；fj干扰信号中心频率
% c=3e8;
Tp=PRI*DutyRatio;
N_Tp=fix(Tp*fs);
% N_PRI=(PRI*fs);
t_Tp=(0:N_Tp-1)/fs;

deltaf=Bj;
filter=fft(fir1(N_Tp-1,deltaf/fs));%滤波器频谱
Noise=ifft(fft(random('Normal',0,sigma,1,N_Tp)).*filter);
Noise_new=zeros(1,N_Tp);
for m=1:N_Tp-1  %Noise作频率，积分得相位Noise_new
    Noise_new(m+1)=Noise(m)+Noise_new(m);
end
FM_jamming_temp=JammingPower*exp(1i*2*pi*(fj*t_Tp+Kfm*Noise_new/fs));
FM_jamming=FM_jamming_temp;
% FM_jamming=zeros(1,N_PRI);
% DelayNumber=ceil(2*R/c*fs);
% FM_jamming(1,DelayNumber:DelayNumber+N_Tp-1)=FM_jamming_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=10e6;K=B/Tp;
f=fj;
J=FM_jamming_temp;
t=(0:length(J)-1)/fs;
s=exp(1i*2*pi*(f*t+1/2*K*t.^2));
h=exp(1i*2*pi*(f*t+1/2*K*t.^2));

s_f=fftshift(fft(s));
h_f=fftshift(fft(h));
j_f=fftshift(fft(J));

s_com_j=fftshift(ifft(j_f.*conj(h_f)));%脉冲压缩信号
s_com=fftshift(ifft(s_f.*conj(h_f)));%真实目标的回波脉压

t1=(0:length(s_com)-1)/fs;
t2=(0:length(s_com_j)-1)/fs;
plot(t1*1e6,abs(real(s_com)),'Linewidth',0.7);
hold on
plot(t2*1e6,abs(real(s_com_j)));
xlabel('时间(us)','fontsize',13,'fontweight','bold');
ylabel('幅度','fontsize',13,'fontweight','bold');
title('脉压','fontsize',13,'fontweight','bold');
legend('目标回波脉压','干扰回波脉压');
end