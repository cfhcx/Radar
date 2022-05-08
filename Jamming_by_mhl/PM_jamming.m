function PM_jamming=PM_jamming(PRI,DutyRatio,JammingPower,sigma,Kpm,Bj,fs,fj,R)
% c=3e8;
Tp=PRI*DutyRatio;
N_Tp=fix(Tp*fs);
% N_PRI=fix(PRI*fs);
t_Tp=(0:N_Tp-1)/fs;

deltaf=Bj;%滤波器截止频率
fil=fft(fir1(N_Tp-1,deltaf/fs));%滤波器频谱
Noise=ifft(fft(random('Normal',0,sigma,1,N_Tp)).*fil);%噪声
% f=(0:length(Noise)-1)/length(Noise)*fs;
% plot(f*1e-6,abs(fft(Noise)));

PM_jamming_temp=JammingPower*exp(1i*2*pi*(fj*t_Tp+Kpm*Noise));
PM_jamming=PM_jamming_temp;
% PM_jamming=zeros(1,N_PRI);
% DelayNumber=ceil(2*R/c*fs);
% PM_jamming(1,DelayNumber:DelayNumber+N_Tp-1)=PM_jamming_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=10e6;K=B/Tp;
f=fj;
J=PM_jamming;
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