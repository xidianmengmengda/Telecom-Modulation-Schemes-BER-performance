%% Part 4

close all; clear;clc
load('A2BPart3.mat')
fc = 20E6;
Ns = 1024;
Rs = 1.5e6;
fs = Rs*Ns;

Nb = 1000;

%split x vector into in phase and quadrature signals


x_qam  = randsrc(1,Nb);
%add value to make even, to be removed at the end
QAM_i = zeros(1,length(x_qam)/4);
QAM_q = QAM_i;


%Assign the 2 bits to a single QAM amplitude
% 00 = -2
% 01 = -1
% 11 = 1
% 10 = 2
% Do for both i and q streams
for i = 1:4:length(x_qam)
   if(x_qam(i) == 1)
       if(x_qam(i+1) == 1)
          QAM_i((i+3)/4) = .33; 
       end
       if(x_qam(i+1) == -1)
           QAM_i((i+3)/4) = 1;
       end
   end

   if(x_qam(i) == -1)
       if(x_qam(i+1) == 1)
           QAM_i((i+3)/4) = -.33;
       end
       if(x_qam(i+1) == -1)
           QAM_i((i+3)/4) = -1;
       end
   end

      if(x_qam(i+2) == 1)
       if(x_qam(i+3) == 1)
          QAM_q((i+3)/4) = .33; 
       end
       if(x_qam(i+3) == -1)
           QAM_q((i+3)/4) = 1;
       end
   end

   if(x_qam(i+2) == -1)
       if(x_qam(i+3) == 1)
           QAM_q((i+3)/4) = -.33;
       end
       if(x_qam(i+3) == -1)
           QAM_q((i+3)/4) = -1;
       end
   end
   
end

t_qam = linspace(0,1/Rs * ceil(Nb/4),Ns*ceil(Nb/4)+1);
t_qam = t_qam(1:end-1);

c_i_QAM = cos(2*pi*fc*t_qam);
c_q_QAM = sin(2*pi*fc*t_qam);

y_i_QAM = c_i_QAM .* kron(QAM_i,ones(1,Ns));
y_q_QAM = c_q_QAM .* kron(QAM_q,ones(1,Ns));


figure();
plot(t_qam,y_i_QAM);

y_QAM = y_i_QAM - y_q_QAM;
figure();
plot(t_qam,y_QAM);

%% Receiver
EbN0_vals2 = -5:16; % Eb/N0 Parameter: from -5dB to 16dB

snr2 = zeros(1,length(EbN0_vals2));
error_QAM = zeros(1,length(EbN0_vals2));

BER_sim2 = [];   % Creates an empty variable that will later collect the BER values for each SNR
for index = 1:length(EbN0_vals2) % sweep through each values of EbNo
    snr2(index) = EbN0_vals2(index)+10*log10(4)-10*log10(Ns); % Convert SNR back to dB as awgn takes dB form
    noisy_y_QAM = awgn(y_QAM,snr2(index),'measured'); % addes AWGN to the modulated signal



% Receiver

    ri_QAM = sum(reshape(noisy_y_QAM .* c_i_QAM,Ns,ceil(Nb/4)));
    rq_QAM = sum(reshape(noisy_y_QAM .* -c_q_QAM,Ns,ceil(Nb/4)));
    r_scatter = [ri_QAM;rq_QAM]';
    %scatterplot(r_scatter);

%plot(ri_QAM);
%Midpoint between the lower and upper amplitudes is 750

%Decode using inverse method from transmitter
% ~1000 = 10
% ~500 = 11
% ~-500 = 01
% ~-1000 = 00
    r_QAM = zeros(1,length(x_qam));
    for i = 1:length(ri_QAM)
        ii = 1+((i-1)*4);
        if(ri_QAM(i) > 333)
            r_QAM(ii) = 1;
            r_QAM(ii+1) = -1;
        end

        if(ri_QAM(i) > 0 && ri_QAM(i) < 333)
            r_QAM(ii) = 1;
            r_QAM(ii+1) = 1;
        end

        if(ri_QAM(i) < 0 && ri_QAM(i) > -333)
            r_QAM(ii) = -1;
            r_QAM(ii+1) = 1;
        end

        if(ri_QAM(i) < -333)
            r_QAM(ii) = -1;
            r_QAM(ii+1) = -1;        
        end


        if(rq_QAM(i) > 333)
            r_QAM(ii+2) = 1;
            r_QAM(ii+3) = -1;
        end

        if(rq_QAM(i) > 0 && rq_QAM(i) < 333)
            r_QAM(ii+2) = 1;
            r_QAM(ii+3) = 1;
        end

        if(rq_QAM(i) < 0 && rq_QAM(i) > -333)
            r_QAM(ii+2) = -1;
            r_QAM(ii+3) = 1;
        end

        if(rq_QAM(i) < -333)
            r_QAM(ii+2) = -1;
            r_QAM(ii+3) = -1;        
        end
    end
    


    error_QAM(index) = sum(r_QAM ~= x_qam);
    BER_sim2 = [BER_sim2,error_QAM(index)/Nb]; % bit error calculation 
end
z = qamdemod(y_QAM,16);

BER_theo2 = (1/4)*3/2*erfc(sqrt(4*0.1*(10.^(EbN0_vals2/10))));
EbN0_vals_lin2 = 10.^(EbN0_vals2/10);
%BER_theo2 = 4/(log2(16))*qfunc(sqrt(3*(log2(16)/(16-1))*EbN0_vals_lin2)); % calculates the theoretical BER 

%BER_theo2=(4/4)*(1-1/sqrt(16))*(1/2)*erfc(sqrt(3*4*(10.^(EbN0_vals2/10))/(16-1))/sqrt(2));

% BER_theo2 = berawgn(EbN0_vals2,'qam',16,);
% plots the theoretical and simulated BER performance 
figure(161), semilogy(EbN0_vals2,BER_sim2,'r+')
hold on
semilogy(EbN0_vals2,BER_theo2,'b')
xlabel('E_b/N_0 dB'), ylabel('BER'), legend('Simulated','Theoretical')
title('Bit Error Rate Performance of BPSK modulation')


%% BER (Justin)


%% 4.2 frequency spectrum
Y_QAM = fft(y_QAM); % fourier transform of the modulated signal
Y_QAM = abs(fftshift(Y_QAM)/fs);    % normalised for plotting
f_QAM = linspace(-fs/2,fs/2,length(Y_QAM)+1); % freqency vector
f_QAM = f_QAM(1:end-1);

figure(114)
plot(f_QAM, Y_QAM)  % plots the freqeucny spectrum of 16-QAM
xlabel('frequency (Hz)'), ylabel('magnitude'), title('Frequency spectrum of 16-QAM')
xlim([15e6 25e6])   % limimts the frequency range
