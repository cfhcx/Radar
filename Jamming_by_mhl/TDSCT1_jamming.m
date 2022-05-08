function TDSCR1_jamming=TDSCT1_jamming(PRI,DutyRatio,JammingPower,B,fs,f,delta_t,a0,R)
%时域采样复制转发 time domain sampling copy repeater
c=3e8;
Tp=PRI*DutyRatio;
% N_PRI=ceil(PRI*fs);
% N_Tp=ceil(Tp*fs);
K=B/Tp;
Delay=2*R/c;
R1=0.1e3;%距离100m
N_Delay=ceil(Delay*fs);
N_delay=ceil(2*R1/c*fs);%目标回波与干扰信号之间的延迟 
delta_N=ceil(delta_t*fs);
N1=ceil(Tp/delta_t);
N_Tp=N1*(delta_N);
N_PRI=ceil(PRI*fs);
% N2=ceil(N_PRI/N_Tp);
t_PRI=(0:N_PRI-1)/fs;
t_Tp=(0:N_Tp-1)/fs;
t=(0:delta_N-1)/fs;
% temp=JammingPower*exp(-1i*2*pi*(f*(t-Delay)+1/2*K*(t-Delay).^2));
rect=repmat(rectpuls(t-a0*delta_t/2,a0*delta_t),1,N1);%在脉宽内进行采样的矩形脉冲串
radar=JammingPower*exp(1i*2*pi*(1/2*K*t_Tp.^2));
radar_sample=rect.*radar;

% radar1=zeros(1,N_PRI);
% radar1(1,1:N_Tp)=radar;
% radar_sample1=zeros(1,N_PRI);
% radar_sample1(1,N_Tp+1:2*N_Tp)=radar_sample;
% figure;subplot(2,1,1);
% plot(t_PRI*1e6,real(radar1));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('截获的雷达信号','fontsize',13,'fontweight','bold');
% subplot(2,1,2);plot(t_PRI*1e6,real(JammingPower*radar_sample1));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('发射的干扰信号','fontsize',13,'fontweight','bold');

TDSCR1_jamming=zeros(1,N_PRI);
TDSCR1_jamming(1,1:N_Tp)=radar_sample;
% TDSCR1_jamming(1,1:N2*N_Tp)=repmat(radar_sample,1,N2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=zeros(1,N_PRI);
J(1,N_Delay+N_Tp:N_Delay+2*N_Tp-1)=radar_sample;

s=zeros(1,N_PRI);
s_temp=exp(1i*2*pi*(1/2*K*t_Tp.^2));
s(1,N_Delay+N_delay+N_Tp:N_Delay+N_delay+2*N_Tp-1)=s_temp;%目标回波信号

% figure;subplot(2,1,1);
% plot(t_PRI*1e6,real(s));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('未处理前的目标回波','fontsize',13,'fontweight','bold');
% subplot(2,1,2);plot(t_PRI*1e6,real(J));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('未处理前的干扰信号','fontsize',13,'fontweight','bold');

h1=exp(1i*2*pi*(1/2*K*t_Tp.^2));
h=conj(fliplr(h1)); %时域匹配滤波为发射信号时间反褶再取共轭 参考信号

s_com_j=conv(J,h); %线性调频信号经过匹配滤波器后的输出(时域卷积)
s_com=conv(s,h);

t1=(0:length(s_com)-1)/fs;
t2=(0:length(s_com_j)-1)/fs;
plot(t1*1e6,abs(real(s_com)),'Linewidth',0.5);
hold on
plot(t2*1e6,abs(real(s_com_j)));
xlim([191 196]);
xlabel('时间','fontsize',13,'fontweight','bold');
ylabel('幅度','fontsize',13,'fontweight','bold');
title('脉压处理后','fontsize',13,'fontweight','bold');
legend('目标回波脉压','干扰回波脉压');
end