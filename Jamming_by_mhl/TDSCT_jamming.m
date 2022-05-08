function TDSCR_jamming=TDSCT_jamming(PRI,DutyRatio,JammingPower,B,fs,f,delta_t,R)
%时域采样复制转发 time domain sampling copy repeater
c=3e8;
Tp=PRI*DutyRatio;
N_PRI=fix(PRI*fs);
N_Tp=fix(Tp*fs);
K=B/Tp;
Delay=2*R/c;
PulseNum=fix((PRI-Tp)/delta_t);%复制的次数
N_start=fix(t_start*fs);
delta_N=fix(delta_t*fs);
N=PulseNum*delta_N;%干扰信号长度
t_Tp=(0:N_Tp-1)/fs;
% t_PRI=(0:N_PRI-1)/fs;
% temp=JammingPower*exp(-1i*2*pi*(f*(t-Delay)+1/2*K*(t-Delay).^2));
radar=JammingPower*exp(1i*2*pi*(1/2*K*t_Tp.^2));
temp=radar(1,N_start+1:N_start+delta_N);
TDSCR_jamming=zeros(1,N_PRI);

% radar1=zeros(1,N_PRI);
% radar1(1,1:N_Tp)=radar;
% temp1=zeros(1,N_PRI);
% temp1(1,N_start+1:N_start+delta_N)=temp;
% figure;subplot(2,1,1);
% plot(t_PRI*1e6,real(radar1));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('雷达信号','fontsize',13,'fontweight','bold');
% subplot(2,1,2);plot(t_PRI*1e6,real(temp1));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('采样','fontsize',13,'fontweight','bold');

TDSCR_jamming(1,delta_N+1:delta_N+N)=JammingPower*repmat(temp,1,PulseNum);
% t=(0:delta_N-1)/fs;
% ff=1e2;
% for k=1:PulseNum
%    TDSCR_jamming(1,(k-1)*delta_N+1:k*delta_N)= JammingPower*exp(1i*2*pi*(1/2*K*t.^2)+1i*2*pi*ff*rand(1));
% end

% figure;plot(t_PRI*1e6,real(TDSCR_jamming));
% xlabel('时间(us)','fontsize',13,'fontweight','bold');
% ylabel('幅度','fontsize',13,'fontweight','bold');
% title('发射干扰信号','fontsize',13,'fontweight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1=0.6e3;%距离600m
N_delay=ceil(2*R1/c*fs);%目标回波与干扰信号之间的延迟 
N_Delay=fix(Delay*fs);

J=zeros(1,N_PRI);
J(1,N_Delay+delta_N+1:N_Delay+delta_N+N)=repmat(temp,1,PulseNum);

s=zeros(1,N_PRI);
s_temp=exp(1i*2*pi*(1/2*K*t_Tp.^2));
s(1,N_Delay+N_delay:N_Delay+N_delay+N_Tp-1)=s_temp;%雷达信号

h1=exp(1i*2*pi*(1/2*K*t_Tp.^2));
% h=zeros(1,N_PRI);
h=conj(fliplr(h1)); %时域匹配滤波为发射信号时间反褶再取共轭 参考信号
% h(1,DelayNumber:DelayNumber+N_Tp-1)=h_temp;%参考信号

s_com_j=conv(J,h); %线性调频信号经过匹配滤波器后的输出(时域卷积)
s_com=conv(s,h);

t1=(0:length(s_com)-1)/fs;
t2=(0:length(s_com_j)-1)/fs;
figure;
plot(t1*1e6,abs(real(s_com)),'Linewidth',0.5);
hold on
plot(t2*1e6,abs(real(s_com_j)));
xlabel('时间(us)','fontsize',13,'fontweight','bold');
ylabel('幅度','fontsize',13,'fontweight','bold');
title('脉压处理后','fontsize',13,'fontweight','bold');
legend('目标回波脉压','干扰回波脉压');
end