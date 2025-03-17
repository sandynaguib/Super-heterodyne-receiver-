% ------------Read audio and process the .wav files----------% 
    [First_message, fs] = audioread('Short_BBCArabic2.wav'); 
    [Second_message, fs] = audioread('Short_RussianVoice.wav');
    maxLength = max([length(First_message), length(Second_message)]);

% ------------------Convert stereo to mono-------------------%
First_message = sum(First_message, 2);
Second_message = sum(Second_message, 2);

%--------Padding signals with zeroes so they have equal length------%
First_message = [First_message; zeros(maxLength - length(First_message), 1)];
Second_message = [Second_message; zeros(maxLength - length(Second_message), 1)];

%------------------Plotting in Time Domain------------------% 

% Create a time axis
time = (0:maxLength-1) / fs;

% Plot the First Message in Time Domain 
figure;
plot(time, First_message);
title('Time Domain of Message 1 BBC Arabic');
xlabel('Time (seconds)');
ylabel('Amplitude');  

% Plot the Second Message in Time Domain
figure;
plot(time, Second_message);
title('Time Domain of Message 2 Russian Voice');
xlabel('Time (seconds)');
ylabel('Amplitude');

%------------Plotting each message in frequency domain--------% 

fft_message_1 = fft(First_message);  % Perform FFT on the message
fft_message_2 = fft(Second_message);
% Creating a frequency axis centered around zero by fftshift 
    fft_shifted_1 = fftshift(fft_message_1);
    fft_shifted_2 = fftshift(fft_message_2);
    Frequency_Range = (-maxLength/2:maxLength/2-1) * fs / maxLength; % Frequency axis from -fs/2 to fs/2  
    % Plot the shifted spectrum (absolute value)
    figure ; 
    plot(Frequency_Range, abs(fft_shifted_1));
    title(['Spectrum of Message 1 BBC Arabic']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    figure ; 
    plot(Frequency_Range, abs(fft_shifted_2));
    title(['Spectrum of Message 2 Russian Voice']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');

%--------------Up Sampling------------------% 
fcarrier_1 = 100e3; % Base carrier frequency (100 KHz)
delta_f = 50e3;  % Frequency increment (50 KHz)

% Increase sampling frequency by a factor to achieve nyquist 
Fs_new = 10 * fs;  % New sampling frequency
% Updated time vector for the new sampling rate

% Interpolate messages to match the new sampling frequency
First_message = interp(First_message, 10);
Second_message = interp(Second_message, 10);

% Updated lengths after interpolation
maxLength_new = max([length(First_message), length(Second_message)]);
t_new = (0:maxLength_new-1) / Fs_new;

% Store the interpolated messages
messages = {First_message, Second_message};

% ---------------Plotting in Time Domain the Upsampled------------------%

%Plot the First Message (Upsampled)
figure;
plot(t_new, First_message);
title('Upsampled Time Domain of Message 1 BBC Arabic');
xlabel('Time (seconds)');
ylabel('Amplitude'); 
 

%Plot the Second Message (Upsampled)
figure;
plot(t_new, Second_message);
title('Upsampled Time Domain of Message 2 Russian Voice');
xlabel('Time (seconds)');
ylabel('Amplitude');

% -------------Plotting in Frequency Domain the Upsampled Signals-------------%

% Compute FFT of upsampled messages
fft_upsampled_1 = fft(First_message);
fft_upsampled_2 = fft(Second_message);

% Create a frequency axis for the new sampling frequency
Frequency_Range_new = (-maxLength_new/2:maxLength_new/2-1) * Fs_new / maxLength_new;

% Shift the FFT results to center at zero frequency
fft_shifted_upsampled_1 = fftshift(fft_upsampled_1);
fft_shifted_upsampled_2 = fftshift(fft_upsampled_2);


%Plot the frequency spectrum of the first upsampled message
figure;
plot(Frequency_Range_new, abs(fft_shifted_upsampled_1));
title('Frequency Spectrum of Upsampled Message 1 (BBC Arabic)');
xlabel('Frequency (Hz)');
ylabel('Magnitude'); 

% Plot the frequency spectrum of the second upsampled message
figure;
plot(Frequency_Range_new, abs(fft_shifted_upsampled_2));
title('Frequency Spectrum of Upsampled Message 2 (Russian Voice)');
xlabel('Frequency (Hz)');
ylabel('Magnitude'); 


%------------------AM Modulation----------------%

% Initialize the FDM signal
FDM_signal = zeros(maxLength_new, 1);

% Iterate over each message and modulate it using DSB-SC
for n = 0:1
    fc = fcarrier_1 + n * delta_f;  % Carrier frequency for the nth signal
    
    % Generate the discrete-time carrier signal at the new sampling rate
    carrier = cos(2 * pi * fc * t_new);  % Discrete-time carrier signal as row vector
    
    % Ensure both message and carrier are the same length
    modulated_signal = messages{n+1} .* carrier'; % Transpose carrier to align with the message

    % Plot the modulated signal
    figure;
    plot(t_new, modulated_signal);
    title(['Modulated Signal ' num2str(n+1) ' at Carrier Frequency ' num2str(fc) ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    % Add to the multiplexed FDM signal
    FDM_signal = FDM_signal + modulated_signal; 
end
%-----------------Plotting FDM--------------------%
%Frequency spectrum of the FDM signal
N = length(FDM_signal);  % Number of samples in the FDM signal
f = (-N/2:N/2-1) * Fs_new / N;  % Frequency axis from -Fs_new/2 to Fs_new/2
spectrum = fftshift(abs(fft(FDM_signal)));  % Compute and shift FFT of the FDM signal

 % Plot the shifted spectrum
figure;
plot(f, spectrum);  
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of FDM Signal');
grid on;

%Time Domain of the FDM signal 
figure;
t_FDM= (0:N-1) /f ; 
plot(t_FDM, FDM_signal);
title('FDM Signal in time domain');
xlabel('Time (seconds)');
ylabel('Amplitude');

% Optionally, you can listen to the FDM signal
soundsc(FDM_signal, Fs_new); 


%-----------------------RF Stage-----------------%

% Design band-pass filter for RF stage
for n = 0:1
    % Center frequency of the filter (carrier frequency)
    fc = fcarrier_1 + n * delta_f;  % Carrier frequency for the nth signal

    % Bandwidth for the filter (choose based on the signal bandwidth)
    BW = 10e3;  % Assume 10 kHz bandwidth for simplicity (adjust as needed)

    % Design the band-pass filter
    BPF = designfilt('bandpassfir','FilterOrder', 100, 'CutoffFrequency1', fc - BW/2, 'CutoffFrequency2', fc + BW/2, ... 
                     'SampleRate', Fs_new);  % Sampling frequency

    % Apply the BPF to the FDM signal to extract the nth signal
    filtered_signal = filter(BPF, FDM_signal);

    % Plot the filtered signal spectrum

    figure;
    spectrum_filtered = fftshift(abs(fft(filtered_signal)));  % Compute FFT and shift
    plot(f, spectrum_filtered);
    title(['Filtered Spectrum at Carrier Frequency ' num2str(fc/1e3) ' kHz']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    %Plot filtered signal in time domain 
    figure ;
    plot(t_fDM,filtered_signal);
    title('Filtered FDM Signal in time domain');
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    

    % Optionally, store the filtered signals for later use
    if n == 0
        RF_signal_1 = filtered_signal;  % First signal
    elseif n == 1
        RF_signal_2 = filtered_signal;  % Second signal
    end
end

%------------------Oscillator and Mixer---------------%

% IF frequency (Intermediate Frequency) = 25 kHz
f_IF = 25e3; 

% Iterate over each filtered signal and apply the oscillator and mixer
for n = 0:1
    % Carrier frequency for the nth signal
    fc = fcarrier_1 + n * delta_f;  
    
    % Frequency of the oscillator is fc + f_IF
    f_oscillator = fc + f_IF;  % Oscillator frequency (sum of fc and IF)
    
    % Generate the oscillator signal (cosine wave with frequency f_oscillator)
    oscillator_signal = cos(2 * pi * f_oscillator * t_new);  % Cosine wave oscillator
    
    % Apply the mixer by multiplying the filtered signal with the oscillator signal
    mixed_signal = RF_signal_1 .* oscillator_signal';  % Mixer: Multiply filtered signal with oscillator
    
    % Store the mixed signal for later stages
    if n == 0
        mixed_signal_1 = mixed_signal;  % Mixed signal for first message
    elseif n == 1
        mixed_signal_2 = mixed_signal;  % Mixed signal for second message
    end
    
    % Plot the frequency domain of the mixed signal
   
    figure;
    mixed_spectrum = fftshift(abs(fft(mixed_signal)));  % Compute and shift FFT of the mixed signal
    plot(f, mixed_spectrum);
    title(['Spectrum of Mixed Signal (Message ' num2str(n+1) ')']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    %Plot the time domain of the mixed signal 
    figure ; 
    plot(t_FDM,mixed_signal); 
    title(['Time Domain of Mixed Signal (Message ' num2str(n+1) ')']);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
   
end
% --------------------IF Stage--------------------%

% Define parameters for the IF band-pass filter
f_center_IF = f_IF;  % Center frequency of the IF stage (25 kHz)
bw_IF = 10e3;  % Bandwidth of the IF filter (adjustable parameter)

% Design the band-pass filter using fdesign.bandpass
IF_filter_design = fdesign.bandpass('N,F3dB1,F3dB2', 10, f_center_IF - bw_IF/2, f_center_IF + bw_IF/2, Fs_new);
IF_filter = design(IF_filter_design, 'butter');

% Apply the IF band-pass filter to the mixed signals
IF_signal_1 = filter(IF_filter, mixed_signal_1);  % Filtered IF signal for message 1
IF_signal_2 = filter(IF_filter, mixed_signal_2);  % Filtered IF signal for message 2

% Compute the spectra of both IF signals
spectrum_IF_signal_1 = fftshift(abs(fft(IF_signal_1)));  % Spectrum of IF signal 1
spectrum_IF_signal_2 = fftshift(abs(fft(IF_signal_2)));  % Spectrum of IF signal 2

% Plot the spectra of both IF signals in separate figures

figure;
subplot(2, 1, 1);
plot(f, spectrum_IF_signal_1);
title('Spectrum of IF Filtered Signal (Message 1)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2, 1, 2);
plot(f, spectrum_IF_signal_2);
title('Spectrum of IF Filtered Signal (Message 2)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% -------------------Baseband Detection-------------------%

% Generate the intermediate frequency carrier
IF_carrier = cos(2 * pi * f_IF * t_new);  % Carrier at IF frequency (25 kHz)

% Mix the IF signals with the IF carrier
baseband_signal_1 = IF_signal_1 .* IF_carrier';  % Transpose to match dimensions
baseband_signal_2 = IF_signal_2 .* IF_carrier';

% Design the low-pass filter (LPF)
f_cutoff = 10e3;  % Cutoff frequency of the LPF (adjustable parameter)
LPF_design = fdesign.lowpass('N,F3dB', 10, f_cutoff, Fs_new);  % Design LPF
LPF = design(LPF_design, 'butter');

% Apply the LPF to obtain baseband signals
detected_signal_1 = filter(LPF, baseband_signal_1);  % Detected baseband signal 1
detected_signal_2 = filter(LPF, baseband_signal_2);  % Detected baseband signal 2

% Plot the detected baseband signals in the time domain

figure;
subplot(2, 1, 1);
plot(t_new, detected_signal_1);
title('Baseband Detected Signal (Message 1)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
plot(t_new, detected_signal_2);
title('Baseband Detected Signal (Message 2)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot the spectra of the detected baseband signals
spectrum_detected_signal_1 = fftshift(abs(fft(detected_signal_1)));
spectrum_detected_signal_2 = fftshift(abs(fft(detected_signal_2)));

figure;
subplot(2, 1, 1);
plot(f, spectrum_detected_signal_1);
title('Spectrum of Baseband Detected Signal (Message 1)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2, 1, 2);
plot(f, spectrum_detected_signal_2);
title('Spectrum of Baseband Detected Signal (Message 2)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% Play the detected baseband signals to verify their correctness
disp('Playing the first detected baseband signal (Message 1)');
soundsc(detected_signal_1, Fs_new); % Play the first detected signal
pause(length(detected_signal_1) / Fs_new + 2); % Pause for playback to complete

disp('Playing the second detected baseband signal (Message 2)');
soundsc(detected_signal_2, Fs_new); % Play the second detected signal
pause(length(detected_signal_2) / Fs_new + 2); % Pause for playback to complete
