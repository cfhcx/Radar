function SMSP_jamming=SMSP_jamming(PRI,DutyRatio,JammingPower,B,fs,f,PulseNum,R)
%PulseNum���������
c=3e8;
Tp=PRI*DutyRatio;
Kj=B/Tp*PulseNum;
N_PRI=PRI*fs;
T_sub=Tp/PulseNum;
N_sub=fix(T_sub*fs); 
N_Tp=N_sub*PulseNum;
t_sub=(0:N_sub-1)/fs;

% SMSP_jamming=zeros(1,N_PRI);
DelayNumber=fix(2*R/c*fs);

% SMSP_jamming_temp1=JammingPower*exp(-1i*2*pi*(fj*t_sub+1/2*Kj*t_sub.^2)); %������
SMSP_jamming_temp1=JammingPower*exp(-1i*2*pi*(1/2*Kj*t_sub.^2)); %������
SMSP_jamming=repmat(SMSP_jamming_temp1,1,PulseNum);
% SMSP_jamming_temp2=repmat(SMSP_jamming_temp1,1,PulseNum);
% SMSP_jamming(1,DelayNumber:DelayNumber+N_Tp-1)=SMSP_jamming_temp2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=SMSP_jamming;
% t_Tp=(0:N_Tp-1)/fs;
t=(0:length(J)-1)/fs;
K=B/Tp;

s=exp(-1i*2*pi*(1/2*K*t.^2));
h=exp(-1i*2*pi*(1/2*K*t.^2));
s_f=fftshift(fft(s));

% h=zeros(1,N_PRI);
% h_temp=s_temp; %ʱ��ƥ���˲�Ϊ�����ź�ʱ�䷴����ȡ���� 
% h(1,1:N_Tp)=h_temp;%�ο��ź�
% h_f=fftshift(fft(h));
h_f=fftshift(fft(h));
j_f=fftshift(fft(J));

s_com_j=fftshift(ifft(j_f.*conj(h_f)));%����ѹ���ź�
s_com=fftshift(ifft(s_f.*conj(h_f)));%��ʵĿ��Ļز���ѹ
% s_com_j=conv(J,h); %���Ե�Ƶ�źž���ƥ���˲���������(ʱ����)
% s_com=conv(s,h);

plot(abs(real(s_com)),'k');
hold on
plot(abs(real(s_com_j)),'r');
xlabel('ʱ��','fontsize',13,'fontweight','bold');
ylabel('����','fontsize',13,'fontweight','bold');
title('���״��źŽ�������ѹ��','fontsize',13,'fontweight','bold');
legend('Ŀ��ز���ѹ','���Żز���ѹ');
end