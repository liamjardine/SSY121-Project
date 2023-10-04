% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency

%message='hello';
%bin= dec2bin(message);
%pack=bin;


function pulse_train = transmitter(pack, fc)

    fs = 26000; %sampling frequency
    R_symb = 150; % Symbol rate
    Q = floor(fs / R_symb); %Q=fsfd, samples per symbol
    fs = Q * R_symb;

    roll_off = 0.35;
    span = 6;

    Tsamp = 1 / fs;
    N = 432; % number of bits to transmit
    assert(length(pack) == N, "The pack to transmitt is of wrong length")
    
    % Constellation or bits to symbol mapping
    const = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 -1i)] / sqrt(2); % Constellation 1 - QPSK/4-QAM

    M = length(const); % Number of symbols in the constellation
    bpsymb = log2(M); % Number of bits per symbol % Symbol rate [symb/s]

    preamble = zadoffChuSeq(859,13)';

    grouped_bits = buffer(pack, bpsymb)'; % Group bits into bits per symbol
    bits_as_symbols = bi2de(grouped_bits, 'left-msb')' + 1; % Bits to symbol index 01=(2^0+2^1)=3(+1 constallation)
    x = const(bits_as_symbols); % Look up symbols using the indices
    x = [preamble x];

    x_upsample = upsample(x, Q);
    % Space the symbols fsfd apart, to enable pulse shaping using conv.

    [pulse, t] = rtrcpuls(roll_off, 1 / R_symb, fs, span);
    pulse_train = conv(pulse, x_upsample);


    tx_signal = pulse_train .* exp(-1i * 2* pi * fc * (0:1:length(pulse_train) - 1) * Tsamp); %carrier
    tx_signal = real(tx_signal);        %take the real part
    tx_signal = tx_signal ./ max(abs(tx_signal)); %normalize it

    player = audioplayer(tx_signal, fs);
    playblocking(player)

end
