function AM_jamming=AM_jamming(PRI,DutyRatio,CarrierPower,NoisePower,Bj,fs,fj,R)
%PRI�����ظ����ڣ�DutyRatioռ�ձȣ�JammingPower���Ź��ʣ�Bj���Ŵ���
%fs����Ƶ�ʣ�fj�����ź�����Ƶ�ʣ�
% c=3e8;
Tp=PRI*DutyRatio;
N_Tp=fix(Tp*fs);
% N_PRI=fix(PRI*fs);
t_Tp=(0:N_Tp-1)/fs;

deltaf=Bj;%�˲�����ֹƵ��
fil=fft(fir1(N_Tp-1,deltaf/fs));%�˲���Ƶ��
Noise=ifft(fft(random('Normal',0,NoisePower,1,N_Tp)).*fil);%���������˲�
AM_jamming_temp=(CarrierPower+Noise).*exp(1i*2*pi*(fj*t_Tp));

% AM_jamming=zeros(1,N_PRI);
% DelayNumber=ceil(2*R/c*fs);
AM_jamming=AM_jamming_temp;
end

