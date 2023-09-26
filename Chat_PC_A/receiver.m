% RECEIVER
%
% This is the receiver structure that you will have to complete.
% The function: receiver(fc) is a setup function for the receiver. Here,
% the audiorecorder object is initialized (see help audiorecorder or
% MATLAB's homepage for more information about the object).
%
% The callback function audioTimerFcn() is a callback function that is
% triggered on a specified time interval (here it is determined in the
% setup function, by the variable time_value)
%
% Your task is to extend this code to work in the project!
%%

function [audio_recorder] = receiver(fc)
    fs = 25000;     %Goal sampling frequency
    R_symb = 400;   %TODO: Choose better wrt frequency mask
    Q = floor(fs/R_symb);   % Samples per symbol
    fs = R_symb * Q;    % Decided sampling frequency, everything is int
    callback_interval = 1; % how often the function should be called in seconds


    assert(fs/2 > fc, "Too low sampling frequency to abide Nyquist.")

    audio_recorder = audiorecorder(fs, 24, 1); % create the recorder

    %attach callback function
    set(audio_recorder, 'TimerPeriod', callback_interval, 'TimerFcn', @audioTimerFcn); % attach a function that should be called every second, the function that is called is specified below.

    %ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
    audio_recorder.UserData.receive_complete = 0; % this is a flag that the while loop in the GUI will check
    audio_recorder.UserData.pack = []; %allocate for data package
    audio_recorder.UserData.pwr_spect = []; %allocate for PSD
    audio_recorder.UserData.const = []; %allocate for constellation
    audio_recorder.UserData.eyed = []; %allocate for eye diagram

    audio_recorder.UserData.fc = fc;
    audio_recorder.UserData.Q = Q;
    audio_recorder.UserData.R_symb = R_symb;



    record(audio_recorder); %start recording
end

% CALLBACK FUNCTION
% This function will be called every [time_value] seconds, where time_value
% is specified above. Note that, as stated in the project MEMO, all the
% fields: pwr_spect, eyed, const and pack need to be assigned if you want
% to get outputs in the GUI.

% So, let's see an example of where we have a pulse train as in Computer
% exercise 2 and let the GUI plot it. Note that we will just create 432
% random bits and hence, the GUI will not be able to decode the message but
% only to display the figures.
% Parameters in the example:
% f_s = 22050 [samples / second]
% R_s = 350 [symbols / second]
% fsfd = f_s/R_s [samples / symbol]
% a = 432 bits
% M = 4 (using QPSK as in computer exercise)

function audioTimerFcn(recObj, event, handles)
    disp('Callback triggered')

    f_sample = recObj.SampleRate;
    Tsample = 1 / f_sample;
    

%%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO: Ansure consistency of roll off with TX
    preamble = [1 2 2 2 3 3 0 0 0]; %TODO: change
    roll_off = 0.3;
    span = 6;
    PA_thresh = 10;     % Placeholder
%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    R_symb = audio_recorder.UserData.R_symb; % Symbol rate [symb/s]    
    T_symb = 1 / R_symb; % Symbol time [s/symb]
    Q = audio_recorder.UserData.Q; % Number of samples per symbol (choose fs such that Q is an integer) [samples/symb]


    f_carrier = recObj.UserData.fc;
    N_bits = 432; % number of bits
    % TODO: ensure consistency of constellation
    const = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 -1i)] / sqrt(2); % Constellation QPSK/4-QAM, [00 01 10 11], GRAY-encoded
    quad_to_bits = [0 0; 1 0; 1 1; 0 1]; % Quadrant number converted to bits

    M = length(const); % Number of symbols in the constellation
    bpsymb = log2(M); % Number of bits per symbol
    N_symbols = N_bits / bpsymb;

    rec_data = getaudiodata(recObj);
    rec_data_downConv = rec_data * exp(-1i * 2 * pi * f_carrier * (0:length(rec_data) - 1) * Tsample);
    rec_data_lowpass = lowpass(rec_data_downConv, f_carrier, f_sample); % Trim LPF if we have noise problems

    preamble_upsample = upsample(const(preamble), Q);
    [pulse,~] = rtrcpuls(roll_off, T_symb, f_sample, span);
    preamble_tx = conv(preamble_upsample, pulse);
    preamble_corr = conv(rec_data_lowpass, fliplr(conj(preamble_tx)));

    [max_correlation, max_index] = max(abs(preamble_corr));
    
    if max_correlation < PA_thresh
        return
    end

    
    data_start_index = max_index + 1;
    data_indices = data_start_index:Q:(data_start_index + (N_symbols - 1) * Q);
    phase_shift = mod(angle(preamble_corr(ind)) * 180 / pi, 360);

    matched_filter = fliplr(conj(pulse));
    %TODO: only convolve on the relevant part of the input stream to speed
    %up
    MF_output = conv(rec_data_lowpass, matched_filter);

    try
        MF_sampled = MF_output(data_indices);
    catch ME
       if ME.identifier == "MATLAB:badsubscript"
           return
       end
    end

    MF_sampled_rotated = MF_sampled .* exp(-1i * (phase_shift / 180) * pi);

    quadrant_number = mod(floor(angle(MF_sampled_rotated) / (pi / 2)), 4) + 1;

    bits = quad_to_bits(quadrant_number, :);
    bits = reshape(bits', 1, []);

    %------------------------------------------------------------------------------
    % HOW TO SAVE DATA FOR THE GUI
    %   NOTE THAT THE EXAMPLE HERE IS ONLY USED TO SHOW HOW TO OUTPUT DATA
    %------------------------------------------------------------------------------

    % Step 1: save the estimated bits
    recObj.UserData.pack = bits;

    % Step 2: save the sampled symbols
    recObj.UserData.const = MF_sampled_rotated;

    % Step 3: provide the matched filter output for the eye diagram (Note here no match filter is used. You should do that)
    recObj.UserData.eyed.r = MF_output;
    recObj.UserData.eyed.fsfd = Q;

    % Step 4: Compute the PSD and save it.
    % !!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
    [pxx, f] = pwelch(rec_data_lowpass, 1024, 768, 1024, f_sample); % note that pwr_spect.f will be normalized frequencies
    f = fftshift(f); %shift to be centered around fs
    f(1:length(f) / 2) = f(1:length(f) / 2) - f_sample; % center to be around zero
    p = fftshift(10 * log10(pxx / max(pxx))); % shift, normalize and convert PSD to dB
    recObj.UserData.pwr_spect.f = f;
    recObj.UserData.pwr_spect.p = p;

    % In order to make the GUI look at the data, we need to set the
    % receive_complete flag equal to 1:
    recObj.UserData.receive_complete = 1;

end
