%% 1) Load FibPho data

signal_465 = [];    % load vector containing 465 nm channel signal values
signal_405 = [];    % load vector containing 405 nm channel signal values
signal_fs = [];     % load FibPho sampling rate (Hz)
first_TTL_FP = [];  % load first TTL timestamp (s)

% remove signal prior to first TTL pulse
onset_FP = first_TTL_FP*signal_fs;

signal_465 = signal_465(onset_FP:end);
signal_405 = signal_405(onset_FP:end);

fs_signal = 1:1:length(signal_465);
sec_signal = fs_signal/signal_fs;   % vector of timestamps (s) for each FP sample point

%% 2) Normalize and plot 

% Determine time interval used for fitting. Use entire trace by defualt and
% check fit. If needed, change start and/or end time of fitting interval.
% This could be in the event of signal artifacts or if signal should only
% be fitted during baseline peridod prior to an intervention

fit_start = 1;                                      % start time (s), must be integer
fit_end = round(sec_signal(end-round(signal_fs)));  % end time (s), must be integer 

fit_interval = (fit_start:1:fit_end); 

MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

% calculate 465 dF/F signal
reg = polyfit(signal_405(round(fit_interval*signal_fs)), signal_465(round(fit_interval*signal_fs)), 1);
a = reg(1);
b = reg(2);
controlFit = a.*signal_405 + b;
normDat = (signal_465 - controlFit)./controlFit;
delta_465 = normDat * 100;  % 465 dF/F signal

% check fit
figure
a = subplot(4,1,1);
plot(sec_signal, signal_405);
title('raw control');
b = subplot(4,1,2);
plot(sec_signal, signal_465);
title('raw signal');
c = subplot(4,1,3);
plot(sec_signal, signal_465);
hold on
plot(sec_signal, controlFit);
title('polynomial fit');
legend({'raw signal','fitted control'})
d = subplot(4,1,4);
plot(sec_signal, delta_465);
title('normalized signal');
xlabel('time (s)');
linkaxes([a,b,c,d],'x');

% filtering trace
delta465_filt = filtfilt(MeanFilter,1,double(delta_465));

% downsampling traces for plotting
ds_factor_FP = 100; % downsampling factor
ds_delta465_filt = downsample(delta465_filt, ds_factor_FP);
ds_sec_signal = downsample(sec_signal, ds_factor_FP);

% Plot filtered trace
figure
plot(ds_sec_signal, ds_delta465_filt)
title('NE (dF/F)');
xlabel('time (s)');

%% 3) loading and plotting EEG and EMG raw data

EEG_rawtrace = [];  % load vector containing EEG channel data
EMG_rawtrace = [];  % load vector containing EMG channeæ data
sampling_freq = []; % load EEG sampling rate (Hz)
first_TTL_EEG = []; % load first TTL timestamp (s)

EEG_time = (1:length(EEG_rawtrace))/sampling_freq; % vector of timestamps (s) for each EEG sample point

% Plot of EEG and EMG traces
figure
a = subplot(2,1,1);
    plot(EEG_time, EMG_rawtrace); 
    ylabel('EMG (V)');
b = subplot(2,1,2);
    plot(EEG_time, EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b],'x');

%% 4) open EEG scoring

%Awake
wake_onset = [];    % onset of each wake bout (s), must be integer
wake_duration = []; % duration of each wake bout (s), must be integer

%Slow-wave sleep
NREM_onset = [];     % onset of each NREM bout (s), must be integer
NREM_duration = [];  % duration of each SWS bout (s), must be integer

%REM
REM_onset = []; % onset of each REM bout (s), must be integer
REM_duration = []; % duration of each REM bout (s), must be integer

% Create binary vectors for each state. Sampling rate is 1 Hz, where idx 1 = time 0-1 s, etc.
wake_binary_vector = zeros(1, round(EEG_time(end))); % vector of zeros matching the length of recording in seconds
for i=1:length(wake_onset)
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1; 
end

NREM_binary_vector = zeros(1, round(EEG_time(end)));
for i=1:length(NREM_onset)
    t = NREM_onset(i)+1; 
    d = NREM_duration(i)-1;
    NREM_binary_vector(t:t+d) = 1;
end

REM_binary_vector =  zeros(1, round(EEG_time(end)));
for i=1:length(REM_onset)
    t = REM_onset(i)+1;
    d = REM_duration(i)-1;
    REM_binary_vector(t:t+d) = 1;
end

% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1;

% check sleep scoring and EEG/EMG data
fig = figure;
a = subplot(2,1,1);
    plot_sleep(EEG_time, EMG_rawtrace, sleepscore_time, wake_binary_vector, NREM_binary_vector, REM_binary_vector);
    ylabel('EMG (V)');
b = subplot(2,1,2);
    plot_sleep(EEG_time, EEG_rawtrace, sleepscore_time, wake_binary_vector, NREM_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b],'x');

% 2-column vectors with on- and offsets for each state
wake_periods = [wake_onset wake_onset+wake_duration];
NREM_periods = [NREM_onset NREM_onset+NREM_duration];
REM_periods = [REM_onset REM_onset+REM_duration];


%% 5) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal (s)
MA_idx = find(wake_duration < MA_maxdur);

% create MA onsets, durations, and binary vectors
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector = zeros(1, round(EEG_time(end)));
for i=1:length(MA_onset)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    MA_binary_vector(t:t+d) = 1;
end

% exclude MAs from wake vectors
wake_woMA_onset = wake_onset;
wake_woMA_onset(MA_idx) = [];
wake_woMA_duration = wake_duration;
wake_woMA_duration(MA_idx) = [];
wake_woMA_binary_vector = zeros(1, round(EEG_time(end)));
for i=1:length(wake_woMA_onset)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

% 2-column vectors with on- and offsets for each state
MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];


%% 6)  Alingment of EEG recording and FP recording

% Cutting EEG/EMG traces leading up to first TTL timestamp
EMG_rawtrace_cut = EMG_rawtrace(round(first_TTL_EEG*sampling_freq):end);
EEG_rawtrace_cut = EEG_rawtrace(round(first_TTL_EEG*sampling_freq):end);
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;

% Cutting part of binary vectors prior to first TTL timestamp
wake_binary_vector_cut = wake_binary_vector(round(first_TTL_EEG+1):end);
NREM_binary_vector_cut = NREM_binary_vector(round(first_TTL_EEG+1):end);
REM_binary_vector_cut = REM_binary_vector(round(first_TTL_EEG+1):end);
MA_binary_vector_cut = MA_binary_vector(round(first_TTL_EEG+1):end);
wake_woMA_binary_vector_cut = wake_woMA_binary_vector(round(first_TTL_EEG+1):end);

% Align onsets, offsets, and durations based on TTL
[wake_onset_cut, wake_offset_cut] = binary_to_OnOff(wake_binary_vector_cut);
wake_duration_cut = wake_offset_cut - wake_onset_cut;

[NREM_onset_cut, NREM_offset_cut] = binary_to_OnOff(NREM_binary_vector_cut);
NREM_duration_cut = NREM_offset_cut - NREM_onset_cut;

[REM_onset_cut, REM_offset_cut] = binary_to_OnOff(REM_binary_vector_cut);
REM_duration_cut = REM_offset_cut - REM_onset_cut;

[MA_onset_cut, MA_offset_cut] = binary_to_OnOff(MA_binary_vector_cut);
MA_duration_cut = MA_offset_cut - MA_onset_cut;

[wake_woMA_onset_cut, wake_woMA_offset_cut] = binary_to_OnOff(wake_woMA_binary_vector_cut);
wake_woMA_duration_cut = wake_woMA_offset_cut - wake_woMA_onset_cut;

% Align period arrays
wake_periods_cut = [wake_onset_cut wake_offset_cut];
NREM_periods_cut = [NREM_onset_cut NREM_offset_cut];
REM_periods_cut = [REM_onset_cut REM_offset_cut];
MA_periods_cut = [MA_onset_cut MA_offset_cut];
wake_woMA_periods_cut = [wake_woMA_onset_cut wake_woMA_offset_cut];

% NB! Here using the cut vector to match the length of binary vectors after being cut to align with FP
sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1;

% 7) Plotting all traces and scorings together
fig = figure;
a = subplot(3,1,1);
    plot_sleep(ds_sec_signal, ds_delta465_filt, sleepscore_time_cut, wake_woMA_binary_vector_cut, NREM_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    title('NE');
b = subplot(3,1,2);
    ds_EEG_time = downsample(EEG_time_cut, 10);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EMG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, NREM_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    ylabel('EMG (V)');
c = subplot(3,1,3);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, NREM_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c],'x');


%% 7) Epoc analysis
% extract mean epoc traces surrounding events of interest (eg. state
% transitions, sensory stimulations, optogenetic manipulations)

epoc_timestamps = MA_onset_cut;   % vector containing timestampts (s) for events of interest (e.g., here micro-arrousal onsets)

% decide length of epoc trace
epoc_start = 30;    % time (s) prior to event
epoc_end = 30;      % time (s) after event

% extract event-related epoc traces
event_465_epoc_collector = [];
for i = 1:length(epoc_timestamps)
    stim_i = epoc_timestamps(i);
    if stim_i > sec_signal(end)-epoc_end % skip if timestamp is too close to end of recording to create full epoc
     continue
    end
    event_465_epoc_i = delta465_filt((stim_i - epoc_start)*signal_fs:(stim_i + epoc_end)*signal_fs);
    event_465_epoc_collector = [event_465_epoc_collector; event_465_epoc_i]; % matrix containing all event-related epoc traces
end

ds_factor = 30;

epoc_FPtime = (1:1:ceil((epoc_start+epoc_end)*signal_fs))/signal_fs-epoc_start;
epoc_FPtime_ds = downsample(epoc_FPtime, ds_factor);

mean_event_465_epocs = downsample(mean(event_465_epoc_collector,1), ds_factor)';

figure
plot(epoc_FPtime_ds, mean_event_465_epocs)
title('mean NE (dF/F)');
xlabel('time from event (s)');

