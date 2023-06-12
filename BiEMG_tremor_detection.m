%%% Function to compute EMG envelope and detect tremor bursts.
% emg_in --> emg raw signal
% fs --> sampling frequency (Hz)
% tremor_th --> tremor detection threshold
% muscle --> muscle to detect tremor

function [EMG_ENVELOPE, BURST_DETECTED] = BiEMG_tremor_detection(emg_in, fs, tremor_th, muscle)

    [BP_b,BP_a] = butter(2, [3 100]/(fs/2), 'bandpass'); % Bandpass filter
    [LP_b,LP_a] = butter(4, 20/(fs/2), 'low'); % Low pass filter
    
    emg_in = emg_in(:,:)'; 
    window_stim = 2*fs;
    
    EMG_fil1 = filtfilt(BP_b,BP_a, emg_in); % BANDPASS filtering
    EMG_fil_h_env = filtfilt(LP_b,LP_a, abs(EMG_fil1)); % ENVELOPE creation
    EMG_fil_h_env = EMG_fil_h_env * 30; % Multipying the increment factor
    tremor_th = tremor_th * max(emg_in(muscle,:));
    for q_muscle=1:muscle
        tremor_burst_loc = [];
        tremor_burst_n = 1;
        for j=2:length(EMG_fil_h_env)-1
            if(EMG_fil_h_env(j,q_muscle) > tremor_th && EMG_fil_h_env(j,q_muscle) < 200 && EMG_fil_h_env(j,q_muscle) > EMG_fil_h_env(j+1,q_muscle) && EMG_fil_h_env(j,q_muscle) > EMG_fil_h_env(j-1,q_muscle))
                tremor_burst_loc(tremor_burst_n) = j;
                tremor_burst_n = tremor_burst_n+1;
            end
        end
        n_bursts = length(tremor_burst_loc);
        ibi = [];
        if(n_bursts>2)
            for i = 1:n_bursts-1
                ibi = [ibi abs(tremor_burst_loc(i) - tremor_burst_loc(i+1))];
            end
            ibi_mean = mean(ibi)/fs;
            freq_tremor = 1/ibi_mean;
            
            next_bursts = [];
            k = 1;
            while((round(tremor_burst_loc(end)/fs + ibi_mean*k*fs) < window_stim))
                next_bursts(k) = round(tremor_burst_loc(end) + ibi_mean*k*fs);
                k = k + 1;
            end
        end
    
        BURST_detected = zeros(length(emg_in),1);
        BURST_detected(tremor_burst_loc) = 100;
        if(isempty(EMG_fil_h_env))
            EMG_fil_h_env = zeros(1000,1);
        end
        BURST_DETECTED(:,q_muscle) = BURST_detected;
    end
    EMG_ENVELOPE = EMG_fil_h_env';
    BURST_DETECTED = BURST_DETECTED';
end