% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency
function transmitter(pack,fc)

fs = 5000; %sampling frequency
tx_signal = randn(1,1e4); %vector of noise
player = audioplayer(tx_signal, fs);       %create an audioplayer object to play the noise at a given sampling frequency
playblocking(player); % Play the noise 
