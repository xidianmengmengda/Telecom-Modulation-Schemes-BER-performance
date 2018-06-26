%% EGB342 Assignment 2B Template 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:/%%%%%%%%%%
%  Run the initialise.m script if you have not already done so.
clear all, close all, clc % clearing and preparing a clean workspace.
%
% Please use variable names identified throughout the brief.
% Where possible, contain code within this script
% i.e. do not create scripts unnecessarily
%
%==========================================================================
% Enter your code below this line:
%==========================================================================
% Student Names: Nick Jenson, Myungjun Choi, Reif Virtue,
% Group Number: 7
%



%% Part 2 - Testing the Wireless Channel
%2.1 INITIALISATION OF THE VARIABLES
fc = 0; % carrier frequency (set as 0 for now - will be modified later)
Rs = 1.5e6; % symbol rate(evaluted in section 2.1)
Ns = 1024; % number of samples/symbol
fs = Ns * Rs;   % sampling rate or pts/Tsym;
Tsym = 1/Rs;    % period of each symbol

% 2.2  LOAD THE MESSAGE DATA ONTO THE WORKSPACE
Nb = 1001; % number of bits in the simulation
load('A2BPart2.mat')    % Load the data into the workspace

% 2.3 CREATE TIME AND FREQUENCY VECTOR
t = linspace(0,Tsym*Nb,Nb*Ns+1);    % the time vector of the message. The vector will start from 
t = t(1:end-1);                     % 0 to Tsym * Ns.

nfft = 2^14;                        % the frequency vector of the message. 
f = linspace(-fs/2,fs/2,nfft+1);    % nfft used to enhance the resolution
f = f(1:end-1);

% 2.4 CREATE A COSINE CARRIER WAVE WITH FREQUENCY OF 20MHZ AND PEAK OF 1V
fc = 20e6   % modify the carrier frequency as required
c = 1*cos(2*pi*fc.*t);  % cosine wave with amplitude of 1 and frequency of 20 MHz


%2.5 MODULATE THE MESSAGE SIGNAL WITH CARRIER WAVE
m = kron(x, ones(1,Ns));    % oversample the message signal so that the vector in lines with the carrier wave 
y= m.*c;    % modulate the message signal on the carrier wave 

figure(1)   % plot the modulated signal
plot(t,y), xlim([0 100e-9]) % display the plot upto first 100ns
xlabel('time (seconds)'), ylabel('voltage (V)'),
title('carrier modulated message signal')

% 2.6 RECOVER THE BIT SEQUENCE USING A CORRELATION RECEIVER 
r = y .* c; % multiplying the basis function with the modulated message signal
r = reshape(r,[Ns,Nb]); % reshape the correlation output so that each row are for each symbol
r = sum(r)/Ns;  % find the mean correlation of each symbol
r = sign(r);    %decision unit - decides whether the symbol is 1 or -1 from its correlation output

% let us compare the recovered bit with the original message visually. 
tbin = t(1:Ns:end); % the time vector for the bits
figure(3)   % plots the recovered and sent bit sequences on the same plot
stem(tbin,r), hold on % recovered bits
stem(tbin,0.5*x,'r') % sent bits (scaled by factor of .5 for better visualisation)
xlabel('Time (sec)'), ylabel('bit')
legend('Recovered (r)','Sent (x)')
title('recovered bit vs original bit') 

biterror = sum(abs(x-r))/2 % bit error calculation

% 2.7 RECOVER THE MESSAGE SIGNAL USING PROVIDED FUNCTION
decoded_message = bits_to_string(r) % converts the bitstream to the message

% 2.8 DETERMINE THE BANDWIDTH OF THE SIGNAL - display the frequency
% spectrum
Y = fft(y,nfft);    % fourier transform of the transmitted message
figure(2);  % plot it with the time vector 
plot(f,abs(fftshift(Y))/nfft);
xlabel('Frequency (Hz)'), ylabel('amplitude'), xlim([0 2e9])
title('BPSK Spectrum for Ns=500 Samples/symbol');



%% 2.12 BIT ERROR PERFORMANCE OF BPSK UNDER AWGN CHANNEL 

EbN0_vals = -5:10; % Eb/N0 Parameter: from -5dB to 10dB
EbN0_vals_lin = 10.^(EbN0_vals/10); % converts the EbNo from dB to linear format
BER_sim = [];    % Creates an empty variable that will later collect the BER values for each SNR

for index = 1:length(EbN0_vals) % sweep through each values of EbNo
    EbN0_lin = EbN0_vals_lin(index); % choose the current value of EbNo 
    snr_lin = EbN0_lin*2*((1/fs)/Tsym); % convets the EbNo to SNR (linear)
    snr = 10*log10(snr_lin); % Convert SNR back to dB as awgn takes dB form
    msg_noisy1 = awgn(y*0.005,snr,'measured'); % addes AWGN to the modulated signal
    
    r2 = msg_noisy1 .* (c); % multiplying the basis function with the modulated message signal (note carrier waveform has been multiplied 
    r = reshape(r2,[Ns,Nb]); % reshape the correlation output so that each row are for each symbol
    r = sum(r)/Ns;  % find the mean correlation of each symbol
    r = sign(r);    %decision unit - decides whether the symbol is 1 or -1 from its correlation output
    BER_sim = [BER_sim,(sum(abs(x-r))/2)/Nb]; % bit error calculation    
end

BER_theo = qfunc(sqrt(2*EbN0_vals_lin)); % calculates the theoretical BER 

% plots the theoretical and simulated BER performance 
figure(5), semilogy(EbN0_vals,BER_sim,'r+') %plots semi-logarithomic scale graph
hold on
semilogy(EbN0_vals,BER_theo,'b')
xlabel('E_b/N_0 dB'), ylabel('BER'), legend('Simulated','Theoretical')
title('Bit Error Rate Performance of BPSK modulation')


