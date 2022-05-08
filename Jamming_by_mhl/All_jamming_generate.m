function jammingsignal=All_jamming_generate(TimeWidth,BandWidth,fs,nrn,nan_num,data_len,mode)
   % TimeWidth 时宽 BandWidth 带宽 fs采样率 nrn距离向采样点数 nan_num脉冲数 data_len干扰个数 mode模式

    if mode==1
    %Jamming 线性函数
        TargetNumber=2;
        TargetVelocity(1:TargetNumber)=[-20 -10];
        TargetDistance(1:TargetNumber)=[3000 4303];
        JammingPower(1:TargetNumber)=[2000 1000];

      %signal=zeros(CPInum,RangeNumber);

        c=3e8;
        Lambda=c/19e6;
        TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda;  
        phase(1:TargetNumber)=2*pi*TargetFd(1:TargetNumber)/nrn;
        DelayNumber(1:TargetNumber)=fix(fs*2*TargetDistance(1:TargetNumber)/c);    
        N=TimeWidth*fs/5;
        f0=0;
        deltaf0=10e6;   
    Ts=1/fs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jammingsignal = zeros(1,data_len);
    for num=1:TargetNumber
        deltaf0=deltaf0+100e6*rand(1);
            k=1-N/2:N-N/2;
                pulse1=JammingPower(num)*exp(j*(pi*(-1.25*BandWidth/TimeWidth)*(k/fs).^2+2*pi*(1*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));
               signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[pulse1];
               jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end
    

%%%%%%%%%%%%%%%%%%%%%============================%%%%%%%%%%%%%%%%%%%%%
elseif mode==2
% Jamming 阶梯波移频
    TargetNumber=4;
    TargetVelocity(1:TargetNumber)=[-20 -10 10 20 ];
    TargetDistance(1:TargetNumber)=[6000 7303 8000 9000];
    JammingPower(1:TargetNumber)=[4000 3000 2000 1000];
  
   

    c=3e8;
    Lambda=c/19e6;
    TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda;  
    phase(1:TargetNumber)=2*pi*TargetFd(1:TargetNumber)/nrn;
    DelayNumber(1:TargetNumber)=fix(fs*2*TargetDistance(1:TargetNumber)/c);    
    N=TimeWidth*fs/5;
    f0=0;
    deltaf0=200e6;   
    Ts=1/fs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jammingsignal = zeros(1,data_len);
    for num=1:2
        deltaf0=deltaf0+100e6*rand(1);
            k=1-N/2:N/4-N/2;
                pulse1=JammingPower(num)*exp(j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(1*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));
            
            k=N/4+1-N/2:N/2-N/2;
                pulse2=JammingPower(num)*exp(j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(2*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));
            
            k=N/2+1-N/2:N/4*3-N/2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                pulse3=JammingPower(num)*exp(j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(3*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));         
            
            k=N/4*3+1-N/2:N-N/2;
                pulse4=JammingPower(num)*exp(j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(4*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));   

               signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[pulse1 pulse2 pulse3 pulse4];
               jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end
    


%%%%%%%%%%%%%%%%%%===============%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif mode==3
%Jamming 分段线性函数
    TargetNumber=2;
    TargetVelocity(1:TargetNumber)=[-20 -10];
    TargetDistance(1:TargetNumber)=[5000 7303];
    JammingPower(1:TargetNumber)=[8000 6000];
  
  %  signal=zeros(CPInum,RangeNumber);

    c=3e8;
    Lambda=c/19e6;
    TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda;  
    phase(1:TargetNumber)=2*pi*TargetFd(1:TargetNumber)/nrn;
    DelayNumber(1:TargetNumber)=fix(fs*2*TargetDistance(1:TargetNumber)/c);    
    N=TimeWidth*fs/5;
    f0=0;
    deltaf0=10e6;   
    Ts=1/fs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jammingsignal = zeros(1,data_len);
    for num=1:2
        deltaf0=deltaf0+100e6*rand(1);
            k=1-N/2:N/4-N/2;
                pulse1=JammingPower(num)*exp(j*(pi*(-1.5*BandWidth/TimeWidth)*(k/fs).^2+2*pi*(1*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));
            
            k=N/4+1-N/2:N/2-N/2;
                pulse2=JammingPower(num)*exp(j*(pi*(-2.5*BandWidth/TimeWidth)*(k/fs).^2+2*pi*(2*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));
            
            k=N/2+1-N/2:N/4*3-N/2;
                pulse3=JammingPower(num)*exp(j*(pi*(-3.5*BandWidth/TimeWidth)*(k/fs).^2+2*pi*(3*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));         
            
            k=N/4*3+1-N/2:N-N/2;
                pulse4=JammingPower(num)*exp(j*(pi*(-0.5*BandWidth/TimeWidth)*(k/fs).^2+2*pi*(4*deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));   

               signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[pulse1 pulse2 pulse3 pulse4];
               jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end
    

elseif mode==4
%Jamming 间歇转发干扰
    TargetNumber=2;
    TargetVelocity(1:TargetNumber)=[-20 -10];
    TargetDistance(1:TargetNumber)=[6000 7103]/2;
    JammingPower(1:TargetNumber)=[4000 6000];% 4000 6000

  %  signal=zeros(CPInum,RangeNumber);

    c=3e8;
    Lambda=c/19e6;
    TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda;  
    phase(1:TargetNumber)=2*pi*TargetFd(1:TargetNumber)./nrn;
    DelayNumber(1:TargetNumber)=fix(fs*2*TargetDistance(1:TargetNumber)/c);    
    N=TimeWidth*fs;
%     deltaf0=10e6;   
    deltaf0=100e6;   
    Ts=TimeWidth;
    tao=0.5*Ts;
    tao_number=fs*tao/5;
    Ts_number=fs*Ts/5;
    temp=zeros(1,N);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jammingsignal = zeros(1,data_len);
    for num=1:TargetNumber
         deltaf0=deltaf0+40e6*rand(1);
            k=1-N/2:N-N/2;
                pulse1=JammingPower(num)*exp(j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(deltaf0+2*TargetFd(num))*(k/fs))+j*(nan_num-1)*phase(num));
            for i=1:5  
                temp(1,(i-1)*fix(Ts_number)+1:(i-1)*fix(Ts_number)+tao_number)=1;
            end
            Signal=pulse1.*temp;
            signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[Signal];
            jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end

elseif mode==5
%Jamming 重复转发干扰
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
    N=TimeWidth*fs/5;
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
    end
end




