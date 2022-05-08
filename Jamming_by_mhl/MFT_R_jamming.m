function MFT_R_jamming=MFT_R_jamming(PRI,DutyRatio,JSR,radarPower,B,fs,f,delta_R,R,n)
%���Ŀ����ƭ֮������ƭ  delta_R�������
JammingPower=10^(JSR/10)*radarPower; 
c=3e8;
Tp=PRI*DutyRatio;
K=B/Tp;
N_Tp=fix(Tp*fs);

delta_t=delta_R*2/c;
delta_N=fix(delta_t*fs);
% t=(0:delta_N-1)/fs;
t_Tp=(0:N_Tp-1)/fs;
N=n*(delta_N+N_Tp);%�����ź��ܳ���

MFT_R_jamming=zeros(1,N);
for k=1:n
    MFT_R_jamming(1,(k-1)*N/n+1:k*N/n)=[JammingPower*exp(1i*2*pi*(f*t_Tp+1/2*K*t_Tp.^2)) zeros(1,delta_N)];
%     MFT_jamming(1,(k-1)*delta_N+1:k*delta_N)=JammingPower.*rectpuls(t-Tp/2,Tp).*exp(-1i*2*pi*((f+fd0+k*delta_Bd)*t+1/2*K*t.^2));
end