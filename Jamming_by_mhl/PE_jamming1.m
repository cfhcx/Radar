function PE_jamming=PE_jamming1(PRI,DutyRatio,JammingPower,B,fs,f0)
T=PRI*DutyRatio;
k=B/T;
Ts=1/fs;
N=ceil(T*fs);
t=linspace(-T/2,T/2,N);
St=exp(1i*2*pi*(f0*t+1/2*k*t.^2)); %�״﷢���ź�  ����Ϊ1
St_noise=awgn(St,10);  %���Ż����յ����ź� �����Ϊ10dB

%%%%%%%%%%%%%%%%%%%%%%�����Ȼ����%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
N1=100;
k_main=[];
f0_main=[];
K=linspace(k-k/10,k+k/10,N1);
F0=linspace(f0-f0/10,f0+f0/10,N1);%���Ʒ�Χ
for time=1:N1
 %��ά����
 for m=1:N1
    for n=1:N1 %�ڵ�Ƶб��Ϊk(m)ʱ����Ƶ�仯
        ML(m,n)=sum(St_noise.*exp(-1i*2*pi*F0(n)*t-1i*pi*K(m)*t.^2));%˲ʱ����ش���
    end
 end
 [C,M]=max(abs(ML(:)));%����ÿһ����ƥ������ǿ��ֵ����λ��
 [A,B]=ind2sub(size(ML),M);%��ǿֵ����λ��()
 k_est=K(A);
 f0_est=F0(B);
 k_main(time)=k_est;
 f0_main(time)=f0_est;
end
k_mean=mean(k_main);
f0_mean=mean(f0_main);
(k-k_mean)/k
(f0-f0_mean)/f0

s=JammingPower*exp(1i*2*pi*(f0_mean*t+1/2*k_mean*t.^2));
PE_jamming=s;

% h1=St;
% h=conj(fliplr(h1)); 
% s_com=conv(s,h);
% St_com=conv(St,h);
% figure;plot(abs(s_com)); hold on; plot(abs(St_com));
% title('��ѹ���');
% xlabel('ʱ���');
% ylabel('����');
% legend('����','��ʵ');

% figure;surf(F0,K,abs(ML));
% title('���ƹ��̿��ӻ�');
% xlabel('Ƶ��');
% ylabel('��Ƶб��');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%�߽�ģ�����㷨%%%%%%%%%%%%%%%%%%%%%%%
if 1
N1=200;
k_main=[];
f0_main=[];
for time=1:N1
 %start to estimate k
 um=[];
 tao=N/2;
 DP2=St_noise.*conj(seqshift(St_noise,t,tao));
 Alpha=linspace(9*tao*pi/(10*T),11*tao*pi/(10*T),N1);
 for m=1:N1
    for n=1:3*N1
        DF2(n)=DP2(n)*exp(-1i*Alpha(m)*t(n));
    end
    um(m)=sum(DF2);
 end
 [m,p]=max(abs(um));
 k_est=Alpha(p)/(2*tao*Ts*pi)
 %end
 %start to estimate f0
 st=St_noise.*exp(-j*k_est*pi*t.^2);
 Beta=linspace(9*f0/10,11*f0/10,N1);
 nm=[];
 for m=1:N1
    for n=1:3*N1
        DF1(n)=st(n)*exp(-j*Beta(m)*t(n));
    end
    nm(m)=sum(DF1);
 end
 [m,r]=max(abs(nm));
 f0_est=Beta(r);
 k_main(time)=k_est;
 f0_main(time)=f0_est;
 %end
end

k_mean=mean(k_main);
f0_mean=mean(f0_main);
(k-k_mean)/k
(f0-f0_mean)/f0

s=JammingPower*exp(1i*2*pi*(f0_mean*t+1/2*k_mean*t.^2));
PE_jamming=s;

h1=St;
h=conj(fliplr(h1)); 
s_com=conv(s,h);
St_com=conv(St,h);
figure;plot(abs(s_com)); hold on; plot(abs(St_com));
% xlim([600 1200]);
title('��ѹ���');
xlabel('ʱ���');
ylabel('����');
legend('����','��ʵ');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function [y,ny]=seqshift(x,nx,k)
    y=x; ny=nx+k;
end