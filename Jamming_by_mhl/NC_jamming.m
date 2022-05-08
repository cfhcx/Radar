function NC_jamming=NC_jamming(PRI,DutyRatio,JammingPower,B,Bj,fs,f,R)
%�������ת�� �����˻�ת������
c=3e8;
Tp=PRI*DutyRatio;
K=B/Tp;
Delay=2*R/c;
N_Delay=fix(Delay*fs);
N_Tp=fix(Tp*fs);
t_Tp=(0:N_Tp-1)/fs;
N_PRI=fix(PRI*fs);

deltaf=Bj;%�˲�����ֹƵ��
fil=fft(fir1(N_Tp-1,deltaf/fs*2));%�˲���Ƶ��
% fil=fftshift(fft(fir1(N_Tp-1,0.1)));%�˲���Ƶ��
Noise_f=fftshift(fft(random('Normal',0,1,1,N_Tp))).*fil;%��˹���������˲������޳��ȸ�˹��������Ƶ��
% Noise_f=fft(random('Normal',0,1,1,N_Tp)).*fil;

% temp1=ifft(Noise_f);
% t=(0:length(temp1)-1)/fs;
% figure;plot(t*1e6,real(temp1));
% xlabel('ʱ��(us)','fontsize',13,'fontweight','bold');
% ylabel('����','fontsize',13,'fontweight','bold');
% title('�����ź�','fontsize',13,'fontweight','bold');

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
% title('Ƶ��ͼ','fontsize',12,'fontweight','bold');
% xlabel('Ƶ��(MHz)','fontsize',12,'fontweight','bold');
% ylabel('��һ������','fontsize',12,'fontweight','bold');
% legend('�����ź�','�״��ź�');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=zeros(1,N_PRI);
s_temp=exp(1i*2*pi*(1/2*K*t_Tp.^2));
s(1,N_Delay:N_Delay+N_Tp-1)=s_temp;%�״��ź�
% s_f=fftshift(fft(s));

% h=zeros(1,N_PRI);
% h_temp=conj( fliplr(s_temp)); %ʱ��ƥ���˲�Ϊ�����ź�ʱ�䷴����ȡ����
% h(1,N_Delay:N_Delay+N_Tp-1)=h_temp;%�ο��ź�
h1=exp(1i*2*pi*(1/2*K*t_Tp.^2));
h=conj(fliplr(h1)); 
% h_f=fftshift(fft(h));
% j_f=fft(NC_jamming);

% s_com_j=ifft(j_f.*h_f);%����ѹ���ź�
% s_com=ifft(s_f.*h_f);%��ʵĿ��Ļز���ѹ
s_com_j=conv(NC_jamming,h); %���Ե�Ƶ�źž���ƥ���˲���������(ʱ����)
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