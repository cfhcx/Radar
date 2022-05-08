function AM_jamming=AM_jamming(PRI,DutyRatio,CarrierPower,NoisePower,Bj,fs,fj,R)
%PRI脉冲重复周期；DutyRatio占空比；JammingPower干扰功率；Bj干扰带宽；
%fs采样频率；fj干扰信号中心频率；
% c=3e8;
Tp=PRI*DutyRatio;
N_Tp=fix(Tp*fs);
% N_PRI=fix(PRI*fs);
t_Tp=(0:N_Tp-1)/fs;

deltaf=Bj;%滤波器截止频率
fil=fft(fir1(N_Tp-1,deltaf/fs));%滤波器频谱
Noise=ifft(fft(random('Normal',0,NoisePower,1,N_Tp)).*fil);%噪声经过滤波
AM_jamming_temp=(CarrierPower+Noise).*exp(1i*2*pi*(fj*t_Tp));

% AM_jamming=zeros(1,N_PRI);
% DelayNumber=ceil(2*R/c*fs);
AM_jamming=AM_jamming_temp;
end

