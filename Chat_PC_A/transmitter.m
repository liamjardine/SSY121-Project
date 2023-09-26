% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency

%message='hello';
%bin= dec2bin(message);
%pack=bin;


function pulse_train = transmitter(pack, fc)

    fs = 20000; %sampling frequency
    R_symb = 400; % Symbol rate
    Q = floor(fs / R_symb); %Q=fsfd, samples per symbol
    fs = Q * R_symb;

    roll_off = 0.35;
    span = 6;

    Tsamp = 1 / fs;
    N = 432; % number of bits to transmit
    % Constellation or bit to symbol mapping
    const = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)] / sqrt(2); % Constellation 1 - QPSK/4-QAM

    M = length(const); % Number of symbols in the constellation
    bpsymb = log2(M); % Number of bits per symbol % Symbol rate [symb/s]
    % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity) [samples/symb]

    %a = randsrc(1, N, [0 1]);
    preamble = [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1];

    a =pack ; % Information bits

    a = [preamble, a];
    m = buffer(a, bpsymb)'; % Group bits into bits per symbol
    m = bi2de(m, 'left-msb')' + 1; % Bits to symbol index 01=(2^0+2^1)=3(+1 constallation)
    x = const(m); % Look up symbols using the indices

    x_upsample = upsample(x, Q);
    % Space the symbols fsfd apart, to enable pulse shaping using conv.

    [pulse, t] = rcpuls(roll_off, 1 / R_symb, fs, span);
    pulse_train = conv(pulse, x_upsample);
    %or maybe fft
%{
npad = length(pulse)+length(x_upsample)-1;
% Calculate the FFT of both signals
pulse=[pulse, zeros(1, npad - length(pulse))];
x_upsample=[x_upsample, zeros(1, npad - length(x_upsample))];
s = ifft(fft(pulse).*fft(x_upsample));
%}

%{
t_vec = Tsamp*(0:1:length(s)-1);
figure; subplot(2,2,1); plot(t_vec,real(s), 'b');
title('real')
subplot(2,2,2); plot(t_vec,imag(s), 'b');
title('imag')
%}

    tx_signal = pulse_train .* exp(-1i * pi * fc * (0:1:length(pulse_train) - 1) * Tsamp); %carrier
    %tx_signal_real=real(s).*sqrt(2)*cos(2*pi*fc*(0:1:length(s)-1)*Tsamp);
    %tx_signal_imag=imag(s).*sqrt(2)*sin(2*pi*fc*(0:1:length(s)-1)*Tsamp);
    %tx_signal=tx_signal_real+tx_signal;
    %tx_signal=real(tx_signal);        %take the real part
    tx_signal = tx_signal ./ max(abs(tx_signal)); %normalize it

    player = audioplayer(tx_signal, fs);
    playblocking(player)

%{
t_vec = Tsamp*(0:1:length(s)-1);
subplot(2,2,3); plot(t_vec+fc,real(tx_signal), 'b');
title('tx_real')
subplot(2,2,4); plot(t_vec+fc,imag(tx_signal), 'b');
title('tx_imag')
%}

end
