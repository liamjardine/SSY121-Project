
function tx_signal = TXdummy(pack, fc)
    fs = 20000; %sampling frequency
    R_symb = 100; % Symbol rate
    Q = floor(fs / R_symb); %Q=fsfd, samples per symbol
    fs = Q * R_symb;
    delay = 0;
    %disp(fs)

    roll_off = 0.35;
    span = 6;

    Tsamp = 1 / fs;
    N = 432; % number of bits to transmit
    %assert(length(pack) == N, "The pack to transmit is of wrong length")
    % Constellation or bit to symbol mapping
    const = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 -1i)] / sqrt(2); % Constellation 1 - QPSK/4-QAM

    M = length(const); % Number of symbols in the constellation
    bpsymb = log2(M); % Number of bits per symbol % Symbol rate [symb/s]
    % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity) [samples/symb]

    %a = randsrc(1, N, [0 1]);
    preamble = [1 1 3 4 4 2 2 2];   % Decimal notation
  
    a = pack ; % Information bits

    %a = [preamble, a];
    m = buffer(a, bpsymb)'; % Group bits into bits per symbol
    m = bi2de(m, 'left-msb')' + 1; % Bits to symbol index 01=(2^0+2^1)=3(+1 constallation)
    msg = [preamble m];
    x = const(msg); % Look up symbols using the indices

    x_upsample = upsample(x, Q);
    % Space the symbols Q apart, to enable pulse shaping using conv.

    [pulse, t] = rtrcpuls(roll_off, 1 / R_symb, fs, span);
    pulse_train = conv(fliplr(pulse), x_upsample);

    tx_signal = pulse_train .* exp(-1i * 2* pi * fc * (0:1:length(pulse_train) - 1) * Tsamp); %carrier
    tx_signal=real(tx_signal);        %take the real part
    tx_signal = tx_signal ./ max(abs(tx_signal)); %normalize it
    tx_signal = [zeros(1,delay) tx_signal];
end