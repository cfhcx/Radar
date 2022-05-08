function RGP_jamming=RGP_jamming(PRI,DutyRatio,JammingPower,B,fs,f,T1,T2,T3,T4,v_t,v,R)
%range gate pull off t_start��ʼʱ�� T1ͣ��ʱ�� T2����ʱ�� T3����ʱ�� v_t�����ٶ�
c=3e8;
Tp=PRI*DutyRatio;
K=B/Tp;
lambda=c/f;
fd=2*v/lambda;
N_Tp=fix(Tp*fs);
N_PRI=fix(PRI*fs);
t_PRI=(0:N_PRI-1)/fs;
Delay=2*R/c;
rect=rectpuls(t_PRI-Tp/2-Delay,Tp);

a=fix(fix(T1*fs)/N_PRI);%ͣ��
N_T1=a*N_PRI;
t1=(0:N_T1-1)/fs;
m=fix(fix(T2*fs)/N_PRI);%����
N_T2=m*N_PRI;
t2=(0:N_T2-1)./fs;
n=fix(fix(T3*fs)/N_PRI);%����
N_T3=n*N_PRI;
t3=(0:N_T3-1)/fs;
h=fix(fix(T4*fs)/N_PRI);%�ر�
N_T4=h*N_PRI;

temp1=repmat(rect.*exp(1i*2*pi*((fd)*(t_PRI-Delay)+1/2*K*(t_PRI-Delay).^2)),1,a);
temp=zeros(N_PRI,m);
for k=0:m-1
    t=t2(k*N_PRI+1:(k+1)*N_PRI);
    s=rectpuls(t_PRI-2*v_t*t/c-Tp/2-Delay,Tp).*exp(1j*2*pi*((fd)*(t_PRI-2*v_t*t/c-Delay)+1/2*K*(t_PRI-2*v_t*t/c-Delay).^2));
    temp(:,k+1)=s;
end
temp2=reshape(temp,1,m*N_PRI);
temp3=repmat(rect.*exp(1i*2*pi*((fd)*(t_PRI-Delay-2*v_t*t2(end)/c)+1/2*K*(t_PRI-2*v_t*t2(end)/c-Delay).^2)),1,n);
temp4=zeros(1,N_T4);

RGP_jamming=JammingPower*[temp1 temp2 temp3 temp4];

echo=repmat(rectpuls(t_PRI-Tp/2-Delay,Tp).*exp(1i*2*pi*((fd)*(t_PRI-Delay)+1/2*K*(t_PRI-Delay).^2)),1,m);
plot(real(echo));
hold on
plot(JammingPower*real(temp2));
title('��ʵĿ���ź�������ź�','fontsize',13,'fontweight','bold');
xlabel('ʱ��','fontsize',13,'fontweight','bold');
ylabel('��ֵ','fontsize',13,'fontweight','bold');
legend('��ʵĿ��','����');
% subplot(2,1,1);
% plot(real(temp2));
% title('�������������ź�ʱ��ͼ','fontsize',12,'fontweight','bold');
% xlabel('ʱ��','fontsize',13,'fontweight','bold');
% ylabel('����','fontsize',13,'fontweight','bold');
% grid on;
% subplot(2,1,2);
% plot(real(echo));
% title('�״�ز��ź�ʱ��ͼ','fontsize',12,'fontweight','bold');
% xlabel('ʱ��','fontsize',13,'fontweight','bold');
% ylabel('����','fontsize',13,'fontweight','bold');
% grid on;

end