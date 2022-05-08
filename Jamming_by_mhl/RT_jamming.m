function RT_jamming=RT_jamming(TimeWidth,BandWidth,fs,nrn,nan_num,data_len)
%重复转发干扰
 TargetNumber=2;
    TargetVelocity(1:TargetNumber)=[-20 -10];
    TargetDistance(1:TargetNumber)=[2000 7303];
    JammingPower(1:TargetNumber)=[5000 7000];

  %  signal=zeros(CPInum,RangeNumber);

    c=3e8;
    Lambda=c/19e6;
    TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda;  
    phase(1:TargetNumber)=2*pi*TargetFd(1:TargetNumber)/nrn;
    DelayNumber(1:TargetNumber)=fix(fs*2*TargetDistance(1:TargetNumber)/c);    
    N=fix(TimeWidth*fs/5);
    deltaf0=10e6;   
   %重复转发参数
    Ts=4e-6;           
    tao=1e-6;
    tao_number=fs*tao/5;
    Ts_number=fs*Ts/5;
    p=zeros(1,N);
    p1=zeros(1,N);
    p2=zeros(1,N);
    p3=zeros(1,N);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jammingsignal = zeros(1,data_len);
    for num=1:TargetNumber
         deltaf0=deltaf0+100e6*rand(1);
            k=1-N/2:N-N/2;
            pulse1=JammingPower(num)*exp(j*(pi*(BandWidth/TimeWidth)*(k/fs).^2+2*pi*(deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));
           %宽带LFM信号
        i=1:N;
        Chirp(1,i)=exp(1j*2*pi*(0.5*BandWidth/TimeWidth*(i*1/fs).^2) );	       %根据公式产生Chirp信号


  
        for jj=1:fix(TimeWidth/Ts)                               %25=TimeWidth/Ts；
            p(1,(jj-1)*fix(Ts_number)+1:(jj-1)*fix(Ts_number)+tao_number)=1;                    
        end

        %重复转发次数=Ts/tao-1,每转发一次增加一个假目标群。本例转发4次，生成三个假目标群。
         p_temp=p.*Chirp;
         p1(1,tao_number+1:N)=p_temp(1,1:fix(N-tao_number));
         p2(1,tao_number+1:N)=p1(1,1:fix(N-tao_number));
         p3(1,tao_number+1:N)=p2(1,1:fix(N-tao_number));
        Signal=p+p1+p2+p3;
            signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[Signal];
            jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end
    RT_jamming=jammingsignal;
end