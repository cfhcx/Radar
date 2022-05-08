function CandI_jamming=CandI_jamming(PRI,DutyRatio,JammingPower,BandWidth,fs,f,RectPulseNum,TimeSlotNum,R)
%RectPulseNum矩形脉冲个数 TimwSlotNum时隙数
c=3e8;
Tp=PRI*DutyRatio;
K=BandWidth/Tp;
N_PRI=fix(PRI*fs);

PulseNum=RectPulseNum*(TimeSlotNum+1);
T_sub=Tp/PulseNum;%子脉宽
N_sub=fix(T_sub*fs);
N_Tp=PulseNum*N_sub;
% t_PRP=0:1/fs:PRP-1/fs;
t_Tp=(0:N_Tp-1)/fs;
% t_sub=(0:N_sub)/fs;

Sample=repmat([ones(1,N_sub),zeros(1,N_sub*TimeSlotNum)],1,RectPulseNum);%矩形取样脉冲
% Radar=exp(-1i*2*pi*(f*t_Tp));
Radar=exp(-1i*2*pi*(1/2*K*t_Tp.^2));
RadarSample=Sample.*Radar;%得到采样后的雷达信号
% figure;plot(t_Tp,Radar);

% CandI_jamming=zeros(1,N_PRI);
CandI_jamming_temp1=zeros(1,N_Tp+TimeSlotNum*N_sub);
for k=0:TimeSlotNum %在时隙中插入子信号
    CandI_jamming_temp1(1,k*N_sub+1:k*N_sub+N_Tp)=CandI_jamming_temp1(1,k*N_sub+1:k*N_sub+N_Tp)+RadarSample;
end
DelayNumber=fix(2*R/c*fs);
CandI_jamming=JammingPower*CandI_jamming_temp1;
% CandI_jamming(1,DelayNumber:DelayNumber+N_Tp-1)=JammingPower*CandI_jamming_temp1(1,1:N_Tp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=CandI_jamming;
t=(0:length(J)-1)/fs;

s=exp(-1i*2*pi*(1/2*K*t.^2));
h=exp(-1i*2*pi*(1/2*K*t.^2));
s_f=fftshift(fft(s));
h_f=fftshift(fft(h));
j_f=fftshift(fft(J));

s_com_j=fftshift(ifft(j_f.*conj(h_f)));%脉冲压缩信号
s_com=fftshift(ifft(s_f.*conj(h_f)));%真实目标的回波脉压

plot(abs(real(s_com)),'Linewidth',1);
hold on
plot(abs(real(s_com_j)));
xlabel('时间','fontsize',13,'fontweight','bold');
ylabel('幅度','fontsize',13,'fontweight','bold');
title('与雷达信号进行脉冲压缩','fontsize',13,'fontweight','bold');
legend('目标回波脉压','干扰回波脉压');
end