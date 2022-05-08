function CandI_jamming=CandI_jamming(PRI,DutyRatio,JammingPower,BandWidth,fs,f,RectPulseNum,TimeSlotNum,R)
%RectPulseNum����������� TimwSlotNumʱ϶��
c=3e8;
Tp=PRI*DutyRatio;
K=BandWidth/Tp;
N_PRI=fix(PRI*fs);

PulseNum=RectPulseNum*(TimeSlotNum+1);
T_sub=Tp/PulseNum;%������
N_sub=fix(T_sub*fs);
N_Tp=PulseNum*N_sub;
% t_PRP=0:1/fs:PRP-1/fs;
t_Tp=(0:N_Tp-1)/fs;
% t_sub=(0:N_sub)/fs;

Sample=repmat([ones(1,N_sub),zeros(1,N_sub*TimeSlotNum)],1,RectPulseNum);%����ȡ������
% Radar=exp(-1i*2*pi*(f*t_Tp));
Radar=exp(-1i*2*pi*(1/2*K*t_Tp.^2));
RadarSample=Sample.*Radar;%�õ���������״��ź�
% figure;plot(t_Tp,Radar);

% CandI_jamming=zeros(1,N_PRI);
CandI_jamming_temp1=zeros(1,N_Tp+TimeSlotNum*N_sub);
for k=0:TimeSlotNum %��ʱ϶�в������ź�
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

s_com_j=fftshift(ifft(j_f.*conj(h_f)));%����ѹ���ź�
s_com=fftshift(ifft(s_f.*conj(h_f)));%��ʵĿ��Ļز���ѹ

plot(abs(real(s_com)),'Linewidth',1);
hold on
plot(abs(real(s_com_j)));
xlabel('ʱ��','fontsize',13,'fontweight','bold');
ylabel('����','fontsize',13,'fontweight','bold');
title('���״��źŽ�������ѹ��','fontsize',13,'fontweight','bold');
legend('Ŀ��ز���ѹ','���Żز���ѹ');
end