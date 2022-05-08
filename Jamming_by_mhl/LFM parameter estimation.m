clc;
clear;
close all;
%%���Ե�Ƶ�ź�
T=10e-6;                                  %p�������ʱ��10us
B=30e6;                                   %���Ե�Ƶ�źŵ�Ƶ�����30MHz
k=B/T;                                      %��Ƶб��
f0=60e6;                                   %��ʼƵ��
Fs=2*B;Ts=1/Fs;                      %����Ƶ�ʺͲ������
N=T/Ts;
t=linspace(-T/2,T/2,N);
St=exp(j*pi*k*t.^2+j*2*pi*f0*t);                    %���Ե�Ƶ�ź�
St_noise=awgn(St,10);                       %�Ӱ׸�˹�������ź�
%%%%%%%%%�����Ȼ����%%%%%%%%%%%%%%%%
% k_main=[];
% f0_main=[];
% for time=1:200
% K=linspace(k-k/10,k+k/10,200);
% F0=linspace(f0-f0/10,f0+f0/10,200);
% for m=1:200
%     for n=1:200
%         ML(m,n)=sum(St_noise.*exp(-j*2*pi*F0(n)*t-j*pi*K(m)*t.^2));
%     end
% end
% [C,I]=max(abs(ML(:)));
% [A,B]=ind2sub(size(ML),I);
% k_est=K(A);
% f0_est=F0(B);
% k_main(time)=k_est;
% f0_main(time)=f0_est;
% end
% k_mean=mean(k_main)
% f0_mean=mean(f0_main)

%%%%%%%%%�߽�ģ���Ⱥ���HAF%%%%%%%%%%%%%%%
k_main=[];
f0_main=[];
for time=1:200
tao=N/2;
DP2=St_noise.*conj(seqshift(St_noise,t,tao));
Alpha=linspace(9*tao*pi/(10*T),11*tao*pi/(10*T),200);
um=[];
for m=1:200
    for n=1:600
        DF2(n)=DP2(n)*exp(-j*Alpha(m)*t(n));
    end
    um(m)=sum(DF2);
end
[m,p]=max(abs(um));
k_est=Alpha(p)/(2*tao*Ts*pi);

st=St_noise.*exp(-j*k_est*pi*t.^2);
Beta=linspace(9*f0/10,11*f0/10,200);
nm=[];
for m=1:200
    for n=1:600
        DF1(n)=st(n)*exp(-j*Beta(m)*t(n));
    end
    nm(m)=sum(DF1);
end
[m,r]=max(abs(nm));
f0_est=Beta(r);
k_main(time)=k_est
f0_main(time)=f0_est
end
k_mean=mean(k_main)
f0_mean=mean(f0_main)

function [y,ny]=seqshift(x,nx,k)
    y=x; ny=nx+k;
end
%%%%%%%%%%%��ʱ����%%%%%%%%%%%%55

%%%%%%%%%%%%%˫��ʱHAF%%%%%%%%%%%











subplot(211)
plot(t*1e6,real(St));
xlabel('ʱ��/us');
title('���Ե�Ƶ�źŵ�ʵ��');
grid on;axis tight;
subplot(212)
freq=linspace(-Fs/2,Fs/2,N);
plot(freq*1e-6,fftshift(abs(fft(St))));
xlabel('Ƶ��/MHz');
title('���Ե�Ƶ�źŵķ�Ƶ����');
grid on;axis tight;

