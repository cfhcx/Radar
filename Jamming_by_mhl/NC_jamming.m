function NC_jamming=NC_jamming(PRI,DutyRatio,JammingPower,B,Bj,fs,f,R)
%噪声卷积转发 噪声乘积转发干扰
c=3e8;
Tp=PRI*DutyRatio;
K=B/Tp;
Delay=2*R/c;
N_Delay=fix(Delay*fs);
N_Tp=fix(Tp*fs);
t_Tp=(0:N_Tp-1)/fs;
N_PRI=fix(PRI*fs);

deltaf=Bj;%滤波器截止频率
fil=fft(fir1(N_Tp-1,deltaf/fs*2));%滤波器频谱
% fil=fftshift(fft(fir1(N_Tp-1,0.1)));%滤波器频谱
Noise_f=fftshift(fft(random('Normal',0,1,1,N_Tp))).*fil;%高斯噪声经过滤波得有限长度高斯白噪声的频谱
% Noise_f=fft(random('Normal',0,1,1,N_Tp)).*fil;

% temp1=ifft(Noise_f);
% t=(0:length(temp1)-1)/fs;
% figure;plot(t*1e6,real(temp1));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('噪声信号','fontsize',13,'fontweight','bold');

radar=exp(1i*2*pi*(1/2*K*t_Tp.^2));
% temp_f=Noise_f.*fft(radar);
% temp=ifft(temp_f);
temp=ifft(Noise_f).*radar;
NC_jamming=zeros(1,N_PRI);
NC_jamming(1,N_Delay:N_Delay+N_Tp-1)=JammingPower*temp;

% radar1=zeros(1,N_PRI);
% radar1(1,1:N_Tp)=radar;
% F1=fftshift(fft(radar1));
% F=fftshift(fft(NC_jamming));
% freq=(0:length(F)-1)*fs/length(F)-B;
% figure;plot(freq*1e-6,abs(F)/max(abs(F)));
% hold on;
% plot(freq*1e-6,abs(F1)/max(abs(F1)));
% title('频谱图','fontsize',12,'fontweight','bold');
% xlabel('频率(MHz)','fontsize',12,'fontweight','bold');
% ylabel('归一化幅度','fontsize',12,'fontweight','bold');
% legend('干扰信号','雷达信号');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=zeros(1,N_PRI);
s_temp=exp(1i*2*pi*(1/2*K*t_Tp.^2));
s(1,N_Delay:N_Delay+N_Tp-1)=s_temp;%雷达信号
% s_f=fftshift(fft(s));

% h=zeros(1,N_PRI);
% h_temp=conj( fliplr(s_temp)); %时域匹配滤波为发射信号时间反褶再取共轭
% h(1,N_Delay:N_Delay+N_Tp-1)=h_temp;%参考信号
h1=exp(1i*2*pi*(1/2*K*t_Tp.^2));
h=conj(fliplr(h1)); 
% h_f=fftshift(fft(h));
% j_f=fft(NC_jamming);

% s_com_j=ifft(j_f.*h_f);%脉冲压缩信号
% s_com=ifft(s_f.*h_f);%真实目标的回波脉压
s_com_j=conv(NC_jamming,h); %线性调频信号经过匹配滤波器后的输出(时域卷积)
s_com=conv(s,h);

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