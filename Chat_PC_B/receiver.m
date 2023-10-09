% RECEIVER
function [audio_recorder] = receiver(f_carrier)
    fs = 28000;     % Goal sampling frequency
    R_symb = 200;   % Symbol rate
    Q = floor(fs / R_symb);     % Samples per symbol
    fs = R_symb * Q;    % Decided sampling frequency, everything is int
    callback_interval = 0.25;   % how often the function should be called in seconds
    min_eye_res = 40;   % Minimum Q to display in eye diagram. Makes GUI somewhat faster


    assert(fs / 2 > f_carrier, "Too low sampling frequency to abide Nyquist.")

    audio_recorder = audiorecorder(fs, 24, 1);  % Create the recorder

    % Attach callback function
    set(audio_recorder, 'TimerPeriod', callback_interval, 'TimerFcn', @audioTimerFcn);

    % ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
    audio_recorder.UserData.receive_complete = 0;   % This is a flag that the while loop in the GUI will check
    audio_recorder.UserData.pack = [];      % Allocate for data package
    audio_recorder.UserData.pwr_spect = [];     % Allocate for PSD
    audio_recorder.UserData.const = [];     % Allocate for constellation
    audio_recorder.UserData.eyed = [];  % Allocate for eye diagram

    audio_recorder.UserData.f_carrier = f_carrier;
    audio_recorder.UserData.Q = Q;
    audio_recorder.UserData.R_symb = R_symb;

    optimum_reduction = ceil(Q/min_eye_res);    % What to divide Q with
    divs = [];
    for i=1:optimum_reduction
        divs = [divs gcd(Q,i)]; % Get divisors without requiring Symbolic Toolbox
    end
    divs = unique(divs);

    % Divide Q with the largest divisor less or equal to or desired one.
    audio_recorder.UserData.best_Q_divisor = max(divs(divs<=optimum_reduction));

    record(audio_recorder); % Start recording
end

% CALLBACK FUNCTION
function audioTimerFcn(recObj, event, handles)
    %disp('Callback triggered')
    t_start = posixtime(datetime);

    f_sample = recObj.SampleRate;
    Tsample = 1 / f_sample;

    %%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    preamble = zadoffChuSeq(859, 13)';
    roll_off = 0.35;    % Roll off factor for root raised cosine
    span = 6;           % Number of T_s to keep of symbol
    PA_thresh = 1;      % Threshold for when to preamble is above noise
    msg_to_keep = 2;    % Number of messages to keep in buffer
    %%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    R_symb = recObj.UserData.R_symb;    % Symbol rate [symb/s]
    T_symb = 1 / R_symb;    % Symbol time [s/symb]
    Q = recObj.UserData.Q;  % Number of samples per symbol (choose fs such that Q is an integer) [samples/symb]
    f_carrier = recObj.UserData.f_carrier;
    N_bits = 432;   % number of bits
    const = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 -1i)] / sqrt(2);   % Constellation QPSK/4-QAM, [00 01 10 11], GRAY-encoded
    quad_to_bits = [0 0; 1 0; 1 1; 0 1];    % Quadrant number converted to bits

    M = length(const);  % Number of symbols in the constellation
    bits_per_symbol = log2(M);
    N_symbols = N_bits / bits_per_symbol;

    % Keep samples equal to {msg_to_keep} messages in buffer
    samples_per_msg = (N_symbols + length(preamble)) / R_symb * f_sample;
    samples_to_keep = ceil(samples_per_msg * msg_to_keep);


    rec_data = getaudiodata(recObj)';
    try
        rec_data = rec_data(1, (end - samples_to_keep):end);
    catch
        %disp("Buffer is shorter than our desired length.")
    end

    rec_data_downConv = rec_data .* exp(1i * 2 * pi * f_carrier .* (0:length(rec_data) - 1) * Tsample);
    rec_data_lowpass = lowpass(rec_data_downConv, f_carrier, f_sample);     % Trim LPF if we have noise problems

    preamble_upsample = upsample(preamble, Q);
    [pulse, ~] = rtrcpuls(roll_off, T_symb, f_sample, span);
    preamble_tx = fftconv(preamble_upsample, pulse);
    preamble_corr = fftconv(rec_data_lowpass, fliplr(conj(preamble_tx)));
    [max_correlation, max_index] = max(abs(preamble_corr));

    %disp(max_correlation + " is max corr now")

    if max_correlation < PA_thresh
        %disp("No PA, only found noise :(")
        return
    end

    data_start_index = max_index + 1;
    data_indices = data_start_index:Q:(data_start_index + (N_symbols - 1) * Q);
    phase_shift = mod(angle(preamble_corr(max_index)) * 180 / pi, 360);

    matched_filter = fliplr(conj(pulse));
    MF_output = fftconv(rec_data_lowpass, matched_filter);

    % Make sure we have the last symbol, which spills over with Q*span
    % samples from the actual sampling point
    if length(MF_output) > (data_indices(end)+Q*(span+1))
        MF_sampled = MF_output(data_indices);
    else
        return
    end

    MF_sampled_rotated = MF_sampled .* exp(-1i * (phase_shift / 180) * pi);

    quadrant_number = mod(floor(angle(MF_sampled_rotated) / (pi / 2)), 4) + 1;

    bits = quad_to_bits(quadrant_number, :);
    bits = reshape(bits', 1, []);



    %%%%% GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The estimated bits
    recObj.UserData.pack = bits;    

    % The sampled symbols (constellation)
    recObj.UserData.const = MF_sampled_rotated / max(abs(MF_sampled_rotated));  
    
    % Eye diagram
    best_Q_divisor = recObj.UserData.best_Q_divisor;    % Downsample Q for faster view of eyediagram
    MF_to_eye = MF_output(data_indices(1):best_Q_divisor:data_indices(end)) .* exp(-1i * (phase_shift / 180) * pi);

    recObj.UserData.eyed.r = MF_to_eye;
    recObj.UserData.eyed.fsfd = Q/best_Q_divisor;


    % PSD 
    % !!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
    [pxx, f] = pwelch(rec_data_lowpass, 1024, 768, 1024, f_sample);     % Note that pwr_spect.f will be normalized frequencies
    f = fftshift(f);    % Shift to be centered around fs
    f(1:length(f) / 2) = f(1:length(f) / 2) - f_sample;     % Center to be around zero
    p = fftshift(10 * log10(pxx / max(pxx)));   % Shift, normalize and convert PSD to dB
    recObj.UserData.pwr_spect.f = f;
    recObj.UserData.pwr_spect.p = p;

    t_end = posixtime(datetime);
    delta = (t_end - t_start) * 1000;
    disp("Data done! Message processing took " + delta + " milliseconds. Peak PA was " + max_correlation)
    recObj.UserData.receive_complete = 1;
end
