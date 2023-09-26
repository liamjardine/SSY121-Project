clf; close all
pack = randi([0 1],1,432);
fc = 2000;


TX = TXdummy(pack,fc);
TX_A = awgn(TX,-8, 'measured');
TX_AR = TX_A*exp(1i *rand(1));



tic
bits = RXdummy(TX_AR,fc);
toc

disp(sum(abs(bits-pack)))


%plot(abs(TX)-abs(TX_A))



%%
TX_down = TX .* exp(-1i * 2 * pi * fc .* (0:length(TX) - 1) * (1/20000));
TX_lowpass = lowpass(TX_down,fc,20000);
plot(real(TX_lowpass)/max(real(TX_lowpass)),"o-")
hold on
plot(real(TX_PULSE)/max(real(TX_PULSE)))
 %%
 TX_down = RXdummy(TX,fc);
 %%
 TX = TXdummy([0 0],fc);
 
 PA_RX = RXdummy(TX,fc);
%%
 plot(abs(TX)/max(abs(TX))-0.05); hold on; plot(abs(PA_RX)/max(abs(PA_RX))); plot(abs(s)/max(abs(s))+0.01); legend(["TX" "RX_PA" "EX5"])

 %% 
%%
%player = audioplayer(TX, 20000);
%playblocking(player)

TX = [TX zeros(1,10000)];
bits = RXdummy(TX,fc);
disp("Done!")
plot(pack - bits)
