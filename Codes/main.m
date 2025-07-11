clc, clearvars;
noisy_speech = load("noisy_speech.txt");
external_noise = load("external_noise.txt");
clean_speech = load("clean_speech.txt");
original_noise = noisy_speech - clean_speech;
fs = 44.1e3;

P_clean_speech = sum((clean_speech).^2);
P_noise_original = sum((original_noise).^2);

filtered_speech_full = noise_cancellation(noisy_speech, external_noise, clean_speech, 0, "full");
[filtered_speech_partial, clean_speech_notched, noisy_speech_notched, filtered_speech_notched] = noise_cancellation(noisy_speech, external_noise, clean_speech, 1e3, "PARTIAL");

snr_original = 10 * log10(P_clean_speech / P_noise_original);

full_suppression_noise = filtered_speech_full - clean_speech;
P_noise_full = sum((full_suppression_noise).^2);
snr_full = 10 * log10(P_clean_speech / P_noise_full);

P_clean_speech_notched = sum((clean_speech_notched).^2);
notched_noise = filtered_speech_notched - clean_speech_notched;
P_new_noise_notched = sum((notched_noise).^2);
snr_partial = 10 * log10(P_clean_speech_notched / P_new_noise_notched);

fprintf("The original SNR is %.2f dB.\n", snr_original); 
fprintf("The SNR after full suppression is %.2f dB.\n", snr_full);
fprintf("The Gain is %.2f dB (for full suppression).\n\n", snr_full - snr_original);

original_noise_notched = noisy_speech_notched - clean_speech_notched;
P_noise_partial_original = sum((original_noise_notched).^2);
snr_partial_original = 10 * log10(P_clean_speech_notched / P_noise_partial_original);

fprintf("The original non-tonal SNR is %.2f dB.\n", snr_partial_original)
fprintf("The SNR after partial suppression is %.2f dB.\n", snr_partial);
fprintf("The Gain is %.2f dB (for partial suppression).\n", snr_partial - snr_partial_original);


% Plots used 

% f = linspace(-fs/2, fs/2, length(clean_speech));
% fft_clean = abs(fftshift(fft(clean_speech)));
% fft_noisy = abs(fftshift(fft(noisy_speech)));
% fft_full = abs(fftshift(fft(filtered_speech_full)));
% fft_partial = abs(fftshift(fft(filtered_speech_partial)));

% subplot(2, 2, 1);
% plot(f, fft_clean);
% grid on;
% title("Spectrum of Clean Speech.");

% subplot(2, 2, 2);
% plot(f, fft_noisy);
% grid on;
% title("Spectrum of Noisy Speech");

% subplot(2, 2, 3);
% plot(f, fft_full);
% grid on;
% title("Spectrum of Output Signal of Full Suppression Mode");

% subplot(2, 2, 4);
% plot(f, fft_partial);
% grid on;
% title("Spectrum of Output Signal of Partial Suppression Mode");
