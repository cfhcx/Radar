
Tp=10e-6;
fj=5e6;
fs=20e6;
N_Tp=ceil(Tp*fs);
t_Tp=linspace(0,Tp,N_Tp);
%  F=fftshift(fft(exp(-1i*2*pi*(fj*t_Tp))));
F=fft(exp(1i*2*pi*(fj*t_Tp)));
funit=(0:length(F)-1)*fs/length(F);
figure;plot(funit,abs(F));