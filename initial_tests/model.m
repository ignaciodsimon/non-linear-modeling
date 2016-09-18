
[clean_audio, sample_rate] = audioread('tube_output_clean.wav');
dist_audio = audioread('tube_output_dist.wav');
% 
% length_audio = round(sample_rate / 1000 * 300);
% clean_audio = clean_audio(151 : 151 + length_audio);
% dist_audio = dist_audio(142 : 142 + length_audio);
% dist_audio = dist_audio / 3 - 0.15;


length_audio = round(sample_rate / 1000 * 1.01);
clean_audio = clean_audio(151 : 151 + length_audio);
dist_audio = dist_audio(142 : 142 + length_audio);
dist_audio = dist_audio / 4;

length_increase_factor = 10000;

interpolated_clean_x = [0 : length_audio / length_increase_factor : length_audio + 1];
interpolated_clean_y = spline([1 : length_audio+1], clean_audio, interpolated_clean_x);

interpolated_dist_x = [0 : length_audio / (length_increase_factor) : length_audio + 1];
interpolated_dist_y = spline([1 : length_audio+1], dist_audio, interpolated_dist_x);

% Beginning trim
interpolated_clean_x = interpolated_clean_x(44 : length(interpolated_clean_x) - 164);
interpolated_clean_y = interpolated_clean_y(44 : length(interpolated_clean_y) - 164);
interpolated_dist_x = interpolated_dist_x(44 : length(interpolated_dist_x) - 164);
interpolated_dist_y = interpolated_dist_y(90 : length(interpolated_dist_y) - 118);

interpolated_clean_y = interpolated_clean_y - interpolated_clean_y(1);
interpolated_dist_y = interpolated_dist_y - interpolated_dist_y(1);


lookup_table = interpolated_dist_y ./ interpolated_clean_y;

plot(lookup_table)
grid
return














plot(interpolated_clean_x, interpolated_clean_y, 'x-');
hold on
plot(interpolated_dist_x, interpolated_dist_y, 'x-');
% xlim([0 3])
% ylim([-0.1 0.3])
grid

return


































% plot(clean_audio)
% hold on
% plot(dist_audio)
% grid
% return

spect_clean = abs(fft(clean_audio));
spect_dist = abs(fft(dist_audio));

phase_clean = phase(fft(clean_audio));
phase_dist = phase(fft(dist_audio));

% plot(phase_clean)
% hold on
% plot(phase_dist)
% grid
% return

spect_clean = spect_clean(1:round(length(spect_clean)/2));
spect_dist = spect_dist(1:round(length(spect_dist)/2));

[fundamental_level, fundamental_pos] = max(spect_dist);
amount_of_harmonics = 20; % <-- this includes the fundamental
harmonics_level = zeros(amount_of_harmonics,1);
phase_delay = zeros(amount_of_harmonics,1);
for i = 1 : amount_of_harmonics
    harmonics_level(i) = spect_dist(((fundamental_pos - 1) * i) +1) / length(spect_dist);
    phase_delay(i) = phase_dist(((fundamental_pos - 1) * i) +1);
end

input_audio = clean_audio;
% input_audio = audioread('bakunin.wav');

% Process clean signal
processed_signal = zeros(size(input_audio));
for i = 1 : amount_of_harmonics
    if mod(phase_delay, 2)
        current_sign = +1;
    else
        current_sign = -1;
    end
    processed_signal = processed_signal + (current_sign * harmonics_level(i) * input_audio.^i);
end

spect_processed = abs(fft(processed_signal));
spect_processed = spect_processed(1 : round(length(spect_processed)/2));
plot(spect_dist)
hold on
plot(spect_processed)
grid
return


% processed_signal = processed_signal / max(abs(processed_signal));
% audiowrite('out.wav', processed_signal, 48000);
plot(input_audio(1000:1200));
hold on
plot(dist_audio(1000:1200))
plot(processed_signal(1000:1200))
grid
return

% plot(20*log10(spect_clean))
% hold on
plot(spect_dist)
% plot(20*log10(spect_dist))
grid
