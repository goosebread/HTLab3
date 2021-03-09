function PPG_reading(sub_num,startTime,endTime,apply_filter,apply_peak)
close all
%folder_name = [num2str(sub_num) '/'];
file_name = ['PPG_Sub' num2str(sub_num) '.txt'];
fid = fopen([file_name],'r');


% read out the data
read_in = textscan(fid,'%f');
data = read_in{1};
% close the file after reading
fclose(fid);

% timestamp of the NI DAQ recording
NI_timestamp(1) = data(1);
NI_timestamp(2) = data(2);
NI_timestamp(3) = data(3);
NI_timestamp(4) = data(4);
NI_timestamp(5) = data(5);
NI_timestamp(6) = data(6);
% remove the timestamp from data
data(1:6) = [];

% sampling rate
sampl_rate = data(1);
data(1) = [];

% number of samples per sampling interval
nr_samples = data(1);
data(1) = [];
ppg_data = data;


% shorten the signal to have the same signal length in all trials
shorten_sig = 1;
%Shortening of the PPG data
ppg_time = [0:1/sampl_rate:numel(ppg_data)/sampl_rate]';
len = min (numel(ppg_data),numel(ppg_time));
ppg_time = ppg_time(1:len);
ppg_data = ppg_data(1:len);
if(shorten_sig == 1) %consider only the the interval given by Start and End Time
    ppg_data = ppg_data(startTime*sampl_rate+1:endTime*sampl_rate);
    
    if startTime~=0
         ppg_time = ppg_time(startTime*sampl_rate+1:endTime*sampl_rate)-ppg_time(startTime*sampl_rate);
    else
        ppg_time = ppg_time(startTime*sampl_rate+1:endTime*sampl_rate);
    end
end

%Display ppg_ir versus time

figure(1)
subplot(3,1,1)
plot(ppg_time, ppg_data),
title('PPG--IR--Raw')
ylabel('Normalized unit')
xlabel('Time (s)')
grid on


if apply_filter == 1
    %%%%%%%%%%%% Lowpass filtering of PPG signal %%%%%%%%%%%%%
    n_order = 5;
    % corner frequency of the lowpass filter in Hz
    f_c = 20;
    ppg_data_filtered = lowpass_filter(ppg_data, sampl_rate, n_order,f_c);
    
    figure(1)
    subplot(3,1,2)
    plot(ppg_time, ppg_data_filtered), hold on
    title('PPG--IR--Filtered')
    ylabel('Normalized unit')
    xlabel('Time (s)')
    grid on
end

if apply_peak == 1
    %%%%%%%%%%%%%% Find peaks %%%%%%%%%%%%%
    max_pulse = 120;
    [ppg_pks,ppg_pks_loc] = find_peaks(ppg_data_filtered, ppg_time, sampl_rate,max_pulse);
    heart_rate = 60*numel(ppg_pks)/(endTime-startTime);
    figure(1)
    subplot(3,1,3)
    plot(ppg_time, ppg_data_filtered), hold on
    plot(ppg_pks_loc, ppg_pks, 'r*')
    title('PPG--IR--Filtered')
    ylabel('Normalized unit')
    xlabel('Time (s)')
    grid on
    legend(['Heart Rate = ' num2str(heart_rate) ' beats/min'])
    
    [m,s]=plotHRV(ppg_pks_loc);
    fprintf("mean is %f\n",m);
    fprintf("standard deviation is %f\n",s);

end
plotPower(ppg_data_filtered,sampl_rate,endTime-startTime);


end


function data_filtered = lowpass_filter(data, sampl_rate, n_order,f_c);
[b_low,a_low] = butter(n_order,f_c/sampl_rate,'low');
data_filtered = filtfilt(b_low,a_low,data);
end

function [data_pks,data_pks_loc] = find_peaks(data, time, sampl_rate,max_pulse);
[pks,pks_loc] = findpeaks(data,'MINPEAKDISTANCE',sampl_rate/(max_pulse/60));
data_pks = pks;
data_pks_loc = time(pks_loc);
end

function [m,s]=plotHRV(hrarray);
%hrarray in seconds
    rates=zeros(1,length(hrarray)-1);%in bpm
        for i=1:length(rates)
            rates(i)=1/(hrarray(i+1)-hrarray(i))*60;
        end
    xbins=min(rates)-0.05:0.1:max(rates)+0.05;
    figure
    histogram(rates,xbins,'Normalization','probability')
    
    m=mean(rates);
    v=var(rates);
    s=v^(0.5);
end

function plotPower(data,sr,time)
    dataP=fft(data);
     figure
    hold on
    freq = (0:(time*sr/2))/time;

    P1=(dataP.*conj(dataP))/(sr*time);
        
        semilogx(freq,P1(1:sr*time/2+1))
        xlabel('Frequency (Hz)')
        xlim([0 3]);
        ylim([0,12]);
        ylabel('Signal Power')
        title(' Power Spectrum')

end

