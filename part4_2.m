%% 4.2 frequency spectrum
Y_QAM = fft(y_QAM); % fourier transform of the modulated signal
Y_QAM = abs(fftshift(Y_QAM)/fs);    % normalised for plotting
f_QAM = linspace(-fs/2,fs/2,length(Y_QAM)+1); % freqency vector
f_QAM = f_QAM(1:end-1);

figure(114)
plot(f_QAM, Y_QAM)  % plots the freqeucny spectrum of 16-QAM
xlabel('frequency (Hz)'), ylabel('magnitude'), title('Frequency spectrum of 16-QAM')
xlim([15e6 25e6])   % limimts the frequency range
