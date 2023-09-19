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
fc = 3000;
fs = 5000; %sampling frequency
audio_recorder = audiorecorder(fs,24,1);% create the recorder

%attach callback function
time_value = 1; % how often the function should be called in seconds
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); % attach a function that should be called every second, the function that is called is specified below.

%ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
audio_recorder.UserData.receive_complete = 0; % this is a flag that the while loop in the GUI will check
audio_recorder.UserData.pack  = []; %allocate for data package
audio_recorder.UserData.pwr_spect = []; %allocate for PSD
audio_recorder.UserData.const = []; %allocate for constellation
audio_recorder.UserData.eyed  = []; %allocate for eye diagram


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

%-----------------------------------------------------------
% THE CODE BELOW IS BASED ON COMPUTER EX 5 AND EX 6:
%-----------------------------------------------------------
disp('Callback triggered')

fs = 5000;                                              % sampling frequency
N = 432;                                                % number of bits
const = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % Constellation 1 - QPSK/4-QAM
M = length(const);                                      % Number of symbols in the constellation
bpsymb = log2(M);                                       % Number of bits per symbol
Rs = 500;                                                % Symbol rate [symb/s]
Ts = 1/Rs;                                              % Symbol time [s/symb]
fsfd = fs/Rs;                                           % Number of samples per symbol (choose fs such that fsfd is an integer) [samples/symb]
bits = randsrc(1,N,[0 1]);                              % Information bits
m_buffer = buffer(bits, bpsymb)';                       % Group bits into bits per symbol
m = bi2de(m_buffer, 'left-msb')'+1;                     % Bits to symbol index
x = const(m);                                           % Look up symbols using the indices
x = awgn(x,15);                                          % add artificial noise
x_upsample = upsample(x, fsfd);                         % Space the symbols fsfd apart, to enable pulse shaping using conv.
span = 6;                                               % Set span = 6
t_vec = -span*Ts: 1/fs :span*Ts;                        % create time vector for one sinc pulse
pulse = sinc(t_vec/Ts);                                 % create sinc pulse with span = 6
pulse_train = conv(pulse,x_upsample);                   % make pulse train



%------------------------------------------------------------------------------
% HOW TO SAVE DATA FOR THE GUI
%   NOTE THAT THE EXAMPLE HERE IS ONLY USED TO SHOW HOW TO OUTPUT DATA
%------------------------------------------------------------------------------

% Step 1: save the estimated bits
recObj.UserData.pack = bits;

% Step 2: save the sampled symbols
recObj.UserData.const = x;

% Step 3: provide the matched filter output for the eye diagram
recObj.UserData.eyed.r = pulse_train;
recObj.UserData.eyed.fsfd = fsfd;

% Step 4: Compute the PSD and save it. 
% !!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
[pxx, f] = pwelch(pulse_train,1024,768,1024, fs); % note that pwr_spect.f will be normalized frequencies
f = fftshift(f); %shift to be centered around fs
f(1:length(f)/2) = f(1:length(f)/2) - fs; % center to be around zero
p = fftshift(10*log10(pxx/max(pxx))); % shift, normalize and convert PSD to dB
recObj.UserData.pwr_spect.f = f;
recObj.UserData.pwr_spect.p = p;

% In order to make the GUI look at the data, we need to set the
% receive_complete flag equal to 1:
recObj.UserData.receive_complete = 1; 
    
end
