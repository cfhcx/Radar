function SWSF_jamming=SWSF_jamming(TimeWidth,BandWidth,fs,nrn,nan_num,data_len)
%½×ÌÝ²¨ÒÆÆµ¸ÉÈÅ
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
                pulse1=JammingPower(num)*exp(1j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(1*deltaf0+2*TargetFd(num))*(k/fs))+1j*(nan_num-1)*phase(num));
            
            k=N/4+1-N/2:N/2-N/2;
                pulse2=JammingPower(num)*exp(1j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(2*deltaf0+2*TargetFd(num))*(k/fs))+1j*(nan_num-1)*phase(num));
            
            k=N/2+1-N/2:N/4*3-N/2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                pulse3=JammingPower(num)*exp(1j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(3*deltaf0+2*TargetFd(num))*(k/fs))+1j*(nan_num-1)*phase(num));         
            
            k=N/4*3+1-N/2:N-N/2;
                pulse4=JammingPower(num)*exp(1j*(pi*(-BandWidth/TimeWidth)*(k/fs).^2+2*pi*(4*deltaf0+2*TargetFd(num))*(k/fs))+1j*(nan_num-1)*phase(num));   

               signaltemp(1,(DelayNumber(num)+1:(DelayNumber(num))+N))=[pulse1 pulse2 pulse3 pulse4];
               jammingsignal=jammingsignal+[signaltemp zeros(1,data_len-length(signaltemp))];
    end
    SWSF_jamming=jammingsignal;
end