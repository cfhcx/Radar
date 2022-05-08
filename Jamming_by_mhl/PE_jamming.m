function PE_jamming=PE_jamming(PRI,DutyRatio,JammingPower,B,fs,f0)
T=PRI*DutyRatio;
k=B/T;
Ts=1/fs;
N=ceil(T*fs);
t=linspace(-T/2,T/2,N);
St=exp(1i*2*pi*(f0*t+1/2*k*t.^2)); %雷达发射信号  幅度为1
% St_noise=awgn(St,10);  %干扰机接收到的信号 信噪比为10dB

%%%%%%%%%%%%%%%%%%%%%%最大似然估计%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
N1=200;
Time=20;
k_main=[];
f0_main=[];
K=linspace(k-k/10,k+k/10,N1);
F0=linspace(f0-f0/10,f0+f0/10,N1);%估计范围
for time=1:Time
 %二维搜索
 for m=1:N1
    for n=1:N1 %在调频斜率为k(m)时，载频变化
        St_noise=awgn(St,10);  %干扰机接收到的信号 信噪比为10dB 共估计Time次，每次St_noise都会因噪声而有些许变化
        ML(m,n)=sum(St_noise.*exp(-1i*2*pi*F0(n)*t-1i*pi*K(m)*t.^2));%瞬时相关处理 exp(-1i*2*pi*F0(n)*t-1i*pi*K(m)*t.^2)作探测信号
    end
 end
 [C,M]=max(abs(ML(:)));%返回每一列中匹配性最强的值及其位置
 [A,B]=ind2sub(size(ML),M);%最强值所在位置()
 k_est=K(A);
 f0_est=F0(B);
 k_main(time)=k_est;
 f0_main(time)=f0_est;
end
k_mean=mean(k_main)
f0_mean=mean(f0_main)
(k-k_mean)/k
(f0-f0_mean)/f0
% figure;plot(1:Time,k_main,'-*');hold on;plot(1:Time,repmat(k_mean,1,Time));
% set(gca,'xtick',1:Time);ylim([4.5e11 5.5e11])
% title('每次的估计值');xlabel('估计序数号');ylabel('估计值');
% figure;plot(1:Time,f0_main,'-*');hold on;plot(1:Time,repmat(f0_mean,1,Time));
% set(gca,'xtick',1:Time);ylim([3.5e7 4.2e7]);title('每次的估计值');xlabel('估计序数号');ylabel('估计值');

s=JammingPower*exp(1i*2*pi*(f0_mean*t+1/2*k_mean*t.^2));
PE_jamming=s;

% h1=St;
% h=conj(fliplr(h1)); 
% s_com=conv(s,h);
% St_com=conv(St,h);
% figure;plot(abs(s_com)); hold on; plot(abs(St_com));
% title('脉压结果');
% xlabel('时间点');
% ylabel('幅度');
% legend('估计','真实');

% figure;surf(F0,K,abs(ML));
% title('估计过程可视化');
% xlabel('频率');
% ylabel('调频斜率');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%高阶模糊度算法%%%%%%%%%%%%%%%%%%%%%%%
if 1
N1=200;
Time=20;
k_main=[];
f0_main=[];
for time=1:Time
 %start to estimate k
 um=[];
 tao=N/2;
 St_noise=awgn(St,10);  %干扰机接收到的信号 信噪比为10dB 
 Alpha=linspace(9*tao*pi/(10*T),11*tao*pi/(10*T),200);
 for m=1:N1
    for n=1:3*N1
        DP2=St_noise.*conj(seqshift(St_noise,t,tao));%对应二阶瞬时矩
        DF2(n)=DP2(n)*exp(-1i*Alpha(m)*t(n));
    end
    um(m)=sum(DF2);
 end
 [m,p]=max(abs(um));
 k_est=Alpha(p)/(2*tao*Ts*pi);
 %end
 %start to estimate f0
 st=St_noise.*exp(-1i*k_est*pi*t.^2);
 Beta=linspace(9*f0/10,11*f0/10,N1);
 nm=[];
 for m=1:10
    for n=1:3*N1
        DF1(n)=st(n)*exp(-1i*2*pi*Beta(m)*t(n));
    end
    nm(m)=sum(DF1);
 end
 [m,r]=max(abs(nm));
 f0_est=Beta(r);
 k_main(time)=k_est;
 f0_main(time)=f0_est;
 %end
end
%对多次估计结果取平均
k_mean=mean(k_main)
f0_mean=mean(f0_main)
(k-k_mean)/k
(f0-f0_mean)/f0
% figure;plot(1:Time,k_main,'-*');hold on;plot(1:Time,repmat(k_mean,1,Time));
% set(gca,'xtick',1:Time);
% title('每次的估计值');xlabel('估计序数号');ylabel('估计值');
s=JammingPower*exp(1i*2*pi*(f0_mean*t+1/2*k_mean*t.^2));
PE_jamming=s;

% h1=St;
% h=conj(fliplr(h1)); 
% s_com=conv(s,h);
% St_com=conv(St,h);
% figure;plot(abs(s_com)); hold on; plot(abs(St_com));
% title('脉压结果');
% xlabel('时间点');
% ylabel('幅度');
% legend('估计','真实');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


