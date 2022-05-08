function VGP_jamming=VGP_jamming(PRI,DutyRatio,JammingPower,B,fs,f,T1,T2,T3,T4,v_f,v,R)
%range gate pull off  
%T1停拖时间 T2拖引时间 T3保持时间 v_f拖引速度 v雷达与干扰机的相对运动速度
c=3e8;
Tp=PRI*DutyRatio;
K=B/Tp;
lambda=c/f;
fd=2*v/lambda;
% N_Tp=fix(Tp*fs);
N_PRI=fix(PRI*fs);
t_PRI=(0:N_PRI-1)/fs;
Delay=2*R/c;
rect=rectpuls(t_PRI-Tp/2-Delay,Tp);
N_T1=fix(fix(T1*fs)/N_PRI)*N_PRI;
t1=(0:N_T1-1)/fs;%停拖
N_T2=fix(fix(T2*fs)/N_PRI)*N_PRI;
t2=(0:N_T2-1)/fs;%拖引
N_T3=fix(fix(T3*fs)/N_PRI)*N_PRI;
t3=(0:N_T3-1)/fs;%保持
N_T4=fix(fix(T4*fs)/N_PRI)*N_PRI;%关闭
%未加载频f
temp1=repmat(rect,1,fix(fix(T1*fs)/N_PRI)).*exp(1i*2*pi*((f+fd)*(t1-Delay)+1/2*K*(t1-Delay).^2));
% temp2=repmat(rect,1,fix(fix(T2*fs)/N_PRI)).*exp(-1i*2*pi*((fd+v_f*(t2-Delay)).*(t2-Delay)+1/2*K*(t2-Delay).^2));
m=fix(fix(T2*fs)/N_PRI);
temp=zeros(N_PRI,m);
for k=0:m-1
    t=t2(k*N_PRI+1:(k+1)*N_PRI);
    s=rectpuls(t_PRI-Tp/2-Delay,Tp).*exp(1i*2*pi*((f+fd+v_f*t).*(t_PRI-Delay)+1/2*K*(t_PRI-Delay).^2));
    temp(:,k+1)=s;
end
temp2=reshape(temp,1,m*N_PRI);
temp3=repmat(rect,1,fix(fix(T3*fs)/N_PRI)).*exp(1i*2*pi*((f+fd+v_f*(N_T2-1)/fs).*(t3-Delay)+1/2*K*(t3-Delay).^2));
temp4=zeros(1,N_T4);

VGP_jamming=JammingPower*[temp1 temp2 temp3 temp4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radar=rect.*exp(1i*2*pi*((f+fd)*(t_PRI-Delay)+1/2*K*(t_PRI-Delay).^2));
radar_f=fft(radar);
s1=JammingPower*temp2(1,(1)*N_PRI:(2)*N_PRI-1);
s1_f=fft(s1);
s2=JammingPower*temp2(1,(m-200)*N_PRI:(m-199)*N_PRI-1);
s2_f=fft(s2);
s3=JammingPower*temp2(1,(m-2)*N_PRI:(m-1)*N_PRI-1);
s3_f=fft(s3);

subplot(3,1,1);
plot(abs(radar_f));
hold on
plot(abs(s1_f));
title('捕获期的某一PRI','fontsize',13,'fontweight','bold');
xlabel('频率','fontsize',13,'fontweight','bold');
ylabel('幅值','fontsize',13,'fontweight','bold');
legend('真实目标','干扰');

subplot(3,1,2);
plot(abs(radar_f));
hold on
plot(abs(s2_f));
title('拖引期的第200个PRI','fontsize',13,'fontweight','bold');
xlabel('频率','fontsize',13,'fontweight','bold');
ylabel('幅值','fontsize',13,'fontweight','bold');
legend('真实目标','干扰');

subplot(3,1,3);
plot(abs(radar_f));
hold on
plot(abs(s3_f));
title('拖引期的第399个PRI','fontsize',13,'fontweight','bold');
xlabel('频率','fontsize',13,'fontweight','bold');
ylabel('幅值','fontsize',13,'fontweight','bold');
legend('真实目标','干扰');
end