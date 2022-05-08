function MFT_V_jamming=MFT_V_jamming(PRI,DutyRatio,JSR,radarPower,B,Bd,fd0,fs,f,delta_R,R,n)
%���Ŀ����ƭ֮�ٶ���ƭ fd0��ʼ������Ƶ�� Bd�����յ��ƴ��� delta_R�������
JammingPower=10^(JSR/10)*radarPower; 
c=3e8;
Tp=PRI*DutyRatio;
K=B/Tp;
N_Tp=fix(Tp*fs);

delta_Bd=Bd/n;
delta_t=delta_R*2/c;
delta_N=fix(delta_t*fs);
% t=(0:delta_N-1)/fs;
t_Tp=(0:N_Tp-1)/fs;
N=n*(delta_N+N_Tp);%�����ź��ܳ���

MFT_V_jamming=zeros(1,N);
for k=1:n
    MFT_V_jamming(1,(k-1)*N/n+1:k*N/n)=[JammingPower*exp(1i*2*pi*((f+fd0+k*delta_Bd)*t_Tp+1/2*K*t_Tp.^2)) zeros(1,delta_N)];
%     MFT_jamming(1,(k-1)*delta_N+1:k*delta_N)=JammingPower.*rectpuls(t-Tp/2,Tp).*exp(-1i*2*pi*((f+fd0+k*delta_Bd)*t+1/2*K*t.^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delay=2*R/c;
N_Delay=fix(Delay*fs);
N_PRI=fix(PRI*fs);
s=zeros(1,N_PRI);
s_temp=exp(1i*2*pi*(1/2*K*t_Tp.^2));
s(1,N_Delay:N_Delay+N_Tp-1)=s_temp;%�״��ź�

s1=exp(1i*2*pi*(1/2*K*t_Tp.^2));
h=conj(fliplr(s1)); %ʱ��ƥ���˲�Ϊ�����ź�ʱ�䷴����ȡ���� �ο��ź�

s_com_j=conv(MFT_V_jamming,h); %���Ե�Ƶ�źž���ƥ���˲���������(ʱ����)
s_com=conv(s,h);

t1=(0:length(s_com)-1)/fs;
t2=(0:length(s_com_j)-1)/fs;
plot(t1*1e6,abs(real(s_com)),'Linewidth',0.7);
hold on
plot(t2*1e6,abs(real(s_com_j)));
xlabel('ʱ��(us)','fontsize',13,'fontweight','bold');
ylabel('����','fontsize',13,'fontweight','bold');
title('��ѹ','fontsize',13,'fontweight','bold');
legend('Ŀ��ز���ѹ','���Żز���ѹ');
end