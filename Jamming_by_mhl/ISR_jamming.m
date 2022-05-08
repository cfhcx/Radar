function ISR_jamming=ISR_jamming(TimeWidth,BandWidth,fs,nrn,nan_num,data_len)
%间歇采样转发干扰 interrupted-sampling and repeater
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
    N=fix(TimeWidth*fs);
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
                temp(1,((i-1)*fix(Ts_number)+1):((i-1)*fix(Ts_number)+tao_number))=1;
            end
            Signal=pulse1.*temp;
            signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[Signal];
            jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end
ISR_jamming=jammingsignal;
end