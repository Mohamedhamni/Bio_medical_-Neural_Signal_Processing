%% EE6209_Assignment2_3336_3426 
%% Part 1
fsignal = bandpass(EEG_RAWdata,[5 20],256);
%fdesign.bandpass
t1  = fsignal(:,1)./fs;
u = 1:2:10;
v = 2:2:10;
figure(1)
for i = 1:5
    val = i + 11;
    subplot(5,2,u(i))
    plot(EEG_RAWdata(:,val),t1);
    strTitle = sprintf('Plot of raw.signal channel %d',val);
    title(strTitle);
    xlabel("Time");
    ylabel("Voltage");
    
    
    subplot(5,2,v(i))
    plot(fsignal(:,val),t1);
    strTitle = sprintf('Plot of filt.signal channel %d',val);
    title(strTitle)
    xlabel("Time");
    ylabel("Voltage");
    
   
end    

raw.fftSignal = fft(EEG_RAWdata);
raw.fftSignal = fftshift(raw.fftSignal);
filt.fftSignal = fft(fsignal);
filt.fftSignal = fftshift(filt.fftSignal);
f = fs/2*linspace(-1,1,fs);

u = 1:2:10;
v = 2:2:10;
figure(2)
for i = 1:5
    val = i + 11;
    subplot(5,2,u(i))
    plot(raw.fftSignal(:,val),t1);
    strTitle = sprintf('Plot of raw.signal magnitude FFT channel %d',val);
    title(strTitle)
    xlabel('Frequency (Hz)');
    ylabel('magnitude');
    
    subplot(5,2,v(i))
    plot(filt.fftSignal(:,val),t1);
    str = sprintf('Plot of filt.signal magnitude FFT channel %d',val);
    title(str)
    xlabel('Frequency (Hz)');
    ylabel('magnitude');
   
end 

%% Part 2
bands = [1 80;81 160;161 240];
t2=linspace(1,768,768);
u = 1:2:10;
v = 2:2:10;
EEG.sum1 =(zeros(768,30));
EEG.sum2 =(zeros(768,30));
EEG.sum3 =(zeros(768,30));
for i =1:3
    for j = bands(i,1):bands(i,2)
        EEG_seg = squeeze(EEG_segmented(j,:,:));
        subplot(3,1,i)
        EEG.("sum"+num2str(i))= EEG.("sum"+num2str(i)) + EEG_seg;
        plot(EEG_seg,t2,'Color','blue')
        strTitle = sprintf('Plot of trails of %d to %d and average of those trails',bands(i,1),bands(i,2));
        title(strTitle)
        xlabel('Time');
        ylabel('magnitude');
        hold on; 
    end
    plot(EEG.("sum"+num2str(i)),t2,'Color','red')
    legend({'all trails','average of trails'},'TextColor','green','FontSize',12)
    legend('boxoff')
    hold off;
end

