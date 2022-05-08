function SF_jamming=SF_jamming(TimeWidth,BandWidth,fs,nrn,nan_num,data_len)
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
    f0=0;%ÔØÆµ ÖÐÐÄÆµÂÊ
    deltaf0=10e6;   
    Ts=1/fs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jammingsignal = zeros(1,data_len);
    for num=1:TargetNumber
        deltaf0=deltaf0+100e6*rand(1);
        k=1-N/2:N-N/2;
        pulse1=JammingPower(num)*exp(1j*(pi*(-1.25*BandWidth/TimeWidth)*(k/fs).^2+2*pi*(1*deltaf0+2*TargetFd(num))*(k/fs))+11j*(nan_num-1)*phase(num));
        signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[pulse1];
        jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end
    SF_jamming=jammingsignal;
end