function bits=RXdummy(TX,f_carrier)
    %disp('Callback triggered')
    fs = 20000; %Goal sampling frequency
    R_symb = 100; %TODO: Choose better wrt frequency mask
    Q = floor(fs / R_symb); % Samples per symbol
    fs = R_symb * Q; % Decided sampling frequency, everything is int

    f_sample = fs;
    Tsample = 1 / f_sample;

    %%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO: Ansure consistency of roll off with TX
    preamble = [1 1 3 4 4 2 2 2]; %TODO: change
    roll_off = 0.35;
    span = 6;
    PA_thresh = 0; % Placeholder
    %%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    T_symb = 1 / R_symb; % Symbol time [s/symb]

    N_bits = 432; % number of bits
    % TODO: ensure consistency of constellation
    const = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 -1i)] / sqrt(2); % Constellation QPSK/4-QAM, [00 01 10 11], GRAY-encoded
    quad_to_bits = [0 0; 1 0; 1 1; 0 1]; % Quadrant number converted to bits

    M = length(const); % Number of symbols in the constellation
    bpsymb = log2(M); % Number of bits per symbol
    N_symbols = N_bits / bpsymb;

    rec_data = TX;
    rec_data_downConv = rec_data .* exp(1i * 2 * pi * f_carrier .* (0:length(rec_data) - 1) * Tsample);
    rec_data_lowpass = lowpass(rec_data_downConv, f_carrier, f_sample); % Trim LPF if we have noise problems
    preamble_upsample = upsample(const(preamble), Q);
    [pulse, ~] = rtrcpuls(roll_off, T_symb, f_sample, span);
    preamble_tx = conv(preamble_upsample, fliplr(conj(pulse)));
    preamble_corr = conv(rec_data_lowpass, fliplr(conj(preamble_tx)));


    [max_correlation, max_index] = max(abs(preamble_corr));
    figure(1)
    plot(abs(preamble_corr))
    title("PA corr")
    if max_correlation < PA_thresh
        disp("No PA, only found noise :(")
        return
    end

    data_start_index = max_index + 1;
    data_indices = data_start_index:Q:(data_start_index + (N_symbols - 1) * Q);
    phase_shift = mod(angle(preamble_corr(max_index)) * 180 / pi, 360);

    matched_filter = fliplr(conj(pulse));
    %TODO: only convolve on the relevant part of the input stream to speed
    %up
    MF_output = conv(rec_data_lowpass, matched_filter);

    try
        MF_sampled = MF_output(data_indices);
    catch ME

        if ME.identifier == "MATLAB:badsubscript"
            disp("nothing yet")
            return
        end

    end

    MF_sampled_rotated = MF_sampled .* exp(-1i * (phase_shift / 180) * pi);
    scatterplot(MF_sampled)

    scatterplot(MF_sampled_rotated)

    quadrant_number = mod(floor(angle(MF_sampled_rotated) / (pi / 2)), 4) + 1;

    bits = quad_to_bits(quadrant_number, :);
    bits = reshape(bits', 1, []);
end
