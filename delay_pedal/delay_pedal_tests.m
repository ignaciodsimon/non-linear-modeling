function delay_pedal_tests()

    input_audio = audioread('dist_marshall_di.wav_processed.wav');

    [ir_short, sample_rate] = audioread('short_delay.wav');
    ir_short = ir_short(1:16000);

    [ir_long, sample_rate] = audioread('long_delay.wav');
    ir_long = ir_long(1:52000);

%     subplot(2,1,1)
%     plot(ir_long);
%     grid
%     subplot(2,1,2)
%     plot(ir_short)
%     grid
%     return

    short_delay_audio = conv(input_audio, ir_short);
    long_delay_audio = conv(input_audio, ir_long);
    short_delay_audio = short_delay_audio(1:length(input_audio));
    long_delay_audio = long_delay_audio(1:length(input_audio));

    left_channel = 2.5 * short_delay_audio;
    right_channel = 2.5 * long_delay_audio;
    
    balance = 20;

    left_channel = left_channel * (100 - balance) / 100;
    right_channel = right_channel * (balance) / 100;

    normalization_value = max([max(abs(left_channel)) max(abs(right_channel))]);
    left_channel = left_channel / normalization_value;
    right_channel = right_channel / normalization_value;
    audiowrite('processed_audio.wav', [left_channel, right_channel], sample_rate);


end
