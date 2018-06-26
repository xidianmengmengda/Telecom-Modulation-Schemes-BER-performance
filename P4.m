%% Part 4
clc; clear all;


load('A2BPart3.mat')
fc = 20E6;
Ns = 1024;
Rs = 1.5e6;
fs = Rs*Ns;

Nb = length(x);

%split x vector into in phase and quadrature signals


x_qam  = [x -1 -1 -1]
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
          QAM_i((i+3)/4) = 1; 
       end
       if(x_qam(i+1) == -1)
           QAM_i((i+3)/4) = 3;
       end
   end

   if(x_qam(i) == -1)
       if(x_qam(i+1) == 1)
           QAM_i((i+3)/4) = -1;
       end
       if(x_qam(i+1) == -1)
           QAM_i((i+3)/4) = -3;
       end
   end

      if(x_qam(i+2) == 1)
       if(x_qam(i+3) == 1)
          QAM_q((i+3)/4) = 1; 
       end
       if(x_qam(i+3) == -1)
           QAM_q((i+3)/4) = 3;
       end
   end

   if(x_qam(i+2) == -1)
       if(x_qam(i+3) == 1)
           QAM_q((i+3)/4) = -1;
       end
       if(x_qam(i+3) == -1)
           QAM_q((i+3)/4) = -3;
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

ri_QAM = sum(reshape(y_QAM .* c_i_QAM,Ns,ceil(Nb/4)));
rq_QAM = sum(reshape(y_QAM .* -c_q_QAM,Ns,ceil(Nb/4)));

plot(ri_QAM);
%Midpoint between the lower and upper amplitudes is 750

%Decode using inverse method from transmitter
% ~1000 = 10
% ~500 = 11
% ~-500 = 01
% ~-1000 = 00
r_QAM = zeros(1,length(x_qam));
for i = 1:length(ri_QAM)
    ii = 1+((i-1)*4);
    if(ri_QAM(i) > 750)
        r_QAM(ii) = 1;
        r_QAM(ii+1) = -1;
    end
    
    if(ri_QAM(i) > 0 && ri_QAM(i) < 750)
        r_QAM(ii) = 1;
        r_QAM(ii+1) = 1;
    end
    
    if(ri_QAM(i) < 0 && ri_QAM(i) > -750)
        r_QAM(ii) = -1;
        r_QAM(ii+1) = 1;
    end

    if(ri_QAM(i) < -750)
        r_QAM(ii) = -1;
        r_QAM(ii+1) = -1;        
    end
    
    
    if(rq_QAM(i) > 1000)
        r_QAM(ii+2) = 1;
        r_QAM(ii+3) = -1;
    end
    
    if(rq_QAM(i) > 0 && rq_QAM(i) < 1000)
        r_QAM(ii+2) = 1;
        r_QAM(ii+3) = 1;
    end
    
    if(rq_QAM(i) < 0 && rq_QAM(i) > -1000)
        r_QAM(ii+2) = -1;
        r_QAM(ii+3) = 1;
    end

    if(rq_QAM(i) < -1000)
        r_QAM(ii+2) = -1;
        r_QAM(ii+3) = -1;        
    end
end

r_QAM = r_QAM(1:length(x));

error_QAM = sum(r_QAM ~= x);


figure();
hold on
stem(r_QAM);
stem(x);



%% BER (Justin)
%% Part 4 - M-QAM Performance in AWGN Channels
clear all; clc

load('A2BPart3.mat')
fc = 20E6;
Ns = 1024;
Rs = 1.5e6;
fs = Rs*Ns;

Nb = length(x);

EbN0_vals2 = -5:16; % Eb/N0 Parameter: from -5dB to 10dB
BER_sim2 = []    % Creates an empty variable that will later collect the BER values for each SNR
x_qam  = [x -1 -1 -1];
%add value to make even, to be removed at the end
QAM_i = zeros(1,length(x_qam)/4);
QAM_q = QAM_i;


%Assign the 2 bits to a single QAM amplitude
% 00 = -2
% 01 = -3
% 11 = 1
% 10 = 3
% Do for both i and q streams

S0 = [-3 -3]; S1 = [-3 -1];
S2 = [-3 3]; S3 = [-3 1];
S4 = [-1 -3]; S5 = [-1 -1];
S6 = [-1 3]; S7 = [-1 1];
S8 = [3 -3]; S9 = [3 -1];
S10 = [3 3]; S11 = [3 1];
S12 = [1 -3]; S13 = [1 -1];
S14= [1 3]; S15 = [1 1];

map16QAM = [S0; S1; S2; S3; S4; S5; S6; S7; S8; S9; S10; S11; S12; S13; S14; S15;]; % Constellation

for i = 1:4:length(x_qam)
   if(x_qam(i) == 1)
       if(x_qam(i+1) == 1)
          QAM_i((i+3)/4) = 1; 
       end
       if(x_qam(i+1) == -1)
           QAM_i((i+3)/4) = 3;
       end
   end

   if(x_qam(i) == -1)
       if(x_qam(i+1) == 1)
           QAM_i((i+3)/4) = -1;
       end
       if(x_qam(i+1) == -1)
           QAM_i((i+3)/4) = -3;
       end
   end

      if(x_qam(i+2) == 1)
       if(x_qam(i+3) == 1)
          QAM_q((i+3)/4) = 1; 
       end
       if(x_qam(i+3) == -1)
           QAM_q((i+3)/4) = 3;
       end
   end

   if(x_qam(i+2) == -1)
       if(x_qam(i+3) == 1)
           QAM_q((i+3)/4) = -1;
       end
       if(x_qam(i+3) == -1)
           QAM_q((i+3)/4) = -3;
       end
   end
   
end

t_qam = linspace(0,1/Rs * ceil(Nb/4),Ns*ceil(Nb/4)+1);
t_qam = t_qam(1:end-1);

c_i_QAM = cos(2*pi*fc*t_qam);
c_q_QAM = sin(2*pi*fc*t_qam);

y_i_QAM = c_i_QAM .* kron(QAM_i,ones(1,Ns));
y_q_QAM = c_q_QAM .* kron(QAM_q,ones(1,Ns));

for index = 1:length(EbN0_vals2) % sweep through each values of EbNo
    snr2 = EbN0_vals2(index)+10*log10(4)-10*log10(Ns); % Convert SNR back to dB as awgn takes dB form
    noisy_y_i_QAM = awgn(y_i_QAM,snr2,'measured'); % addes AWGN to the modulated signal
    noisy_y_q_QAM = awgn(y_q_QAM,snr2,'measured'); % addes AWGN to the modulated signal
    
%     figure();
%     plot(t_qam,y_i_QAM);

    y_QAM = noisy_y_i_QAM - noisy_y_q_QAM;
%     figure();
%     plot(t_qam,y_QAM);

% Receiver

    ri_QAM = sum(reshape(y_QAM,Ns,ceil(Nb/4)));
    rq_QAM = sum(reshape(y_QAM,Ns,ceil(Nb/4)));

    % Compute distance metrics
    metrics = NaN*zeros(length(map16QAM),length(x_qam)/4);
    for jj = 1:16
        metrics(jj,:) = (ri_QAM-map16QAM(jj,1)).^2+(rq_QAM-map16QAM(jj,2)).^2; % ML detection
    end
    
    [min_metric, location] =min(metrics) ; % Find the index of the minimum metric
    RxBits = decimalToBinaryVector(location-1,4).'; % Convert the M-ary received data to binary
    RxBits = RxBits(:).';
    RxBits(RxBits == 0) = -1;
    % Count the bit errors
    berf(index) = sum(x_qam~=RxBits);
    
figure(index)    
plot(ri_QAM);
%Midpoint between the lower and upper amplitudes is 750

%Decode using inverse method from transmitter
% ~1000 = 10
% ~500 = 11
% ~-500 = 01
% % ~-1000 = 00
%     r_QAM = zeros(1,length(x_qam));
%     for i = 1:length(ri_QAM)
%         ii = 1+((i-1)*4);
%         if(ri_QAM(i) > 3*max(ri_QAM)/4)
%             r_QAM(ii) = 1;
%             r_QAM(ii+1) = -1;
%         end
% 
%         if(ri_QAM(i) > 0 && ri_QAM(i) < 2*max(ri_QAM)/4)
%             r_QAM(ii) = 1;
%             r_QAM(ii+1) = 1;
%         end
% 
%         if(ri_QAM(i) < 0 && ri_QAM(i) > 1*max(ri_QAM)/4)
%             r_QAM(ii) = -1;
%             r_QAM(ii+1) = 1;
%         end
% 
%         if(ri_QAM(i) < 1*max(ri_QAM)/4)
%             r_QAM(ii) = -1;
%             r_QAM(ii+1) = -1;        
%         end
% 
% 
%         if(rq_QAM(i) > 1000)
%             r_QAM(ii+2) = 1;
%             r_QAM(ii+3) = -1;
%         end
% 
%         if(rq_QAM(i) > 0 && rq_QAM(i) < 1000)
%             r_QAM(ii+2) = 1;
%             r_QAM(ii+3) = 1;
%         end
% 
%         if(rq_QAM(i) < 0 && rq_QAM(i) > -1000)
%             r_QAM(ii+2) = -1;
%             r_QAM(ii+3) = 1;
%         end
% 
%         if(rq_QAM(i) < -1000)
%             r_QAM(ii+2) = -1;
%             r_QAM(ii+3) = -1;        
%         end
%     end

%     r_QAM = r_QAM(1:length(x));
% 
%     error_QAM = sum(r_QAM ~= x);
%     BER_sim2 = [BER_sim2,error_QAM/Nb]; % bit error calculation 
end

BER_sim2 = berf/(length(x_qam));% k bits per symbol - bit error

EbN0_vals_lin2 = 10  .^(EbN0_vals2/10);
BER_theo2 = (log2(16)/2)*qfunc(sqrt(3*(log2(16)/(16-1))*EbN0_vals_lin2)); % calculates the theoretical BER 

% plots the theoretical and simulated BER performance 
figure(15), semilogy(EbN0_vals2,BER_sim2,'r+')
hold on
semilogy(EbN0_vals2,BER_theo2,'b')
xlabel('E_b/N_0 dB'), ylabel('BER'), legend('Simulated','Theoretical')
title('Bit Error Rate Performance of BPSK modulation')

%% frequency

Y_QAM = fft(y_QAM);
Y_QAM = abs(fftshift(Y_QAM)/fs);
f_QAM = linspace(-fs/2,fs/2,length(Y_QAM)+1);
f_QAM = f_QAM(1:end-1);
figure()
plot(f_QAM,Y_QAM);