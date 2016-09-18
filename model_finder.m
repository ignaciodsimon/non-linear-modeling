function model_finder()


    [mls_input, sample_rate] = audioread('tube_amp_io_recordings/tube_amp_mls_in.wav');
    mls_output = audioread('tube_amp_io_recordings/tube_amp_mls_out.wav');
    sine_input = audioread('tube_amp_io_recordings/tube_amp_tone_dist_in.wav');
    sine_output = audioread('tube_amp_io_recordings/tube_amp_tone_dist_out.wav');


    % 1 - Find the linear characteristics: band filtering shape

    % (The MLS recordings are assumed to be MLS signals looped).
    ir_band_filter = xcorr(mls_output, mls_input);
    ir_band_filter = ir_band_filter(round(length(ir_band_filter)/2) : round(length(ir_band_filter)/2) + 2300);
    ir_band_filter = ir_band_filter / max(abs(ir_band_filter));

    % To correct for a possible phase inversion
    [~, max_pos] = max(abs(ir_band_filter));
    if ir_band_filter(max_pos) < 0
        ir_band_filter = -ir_band_filter;
    end


%     disp('Convolving ...')
%     input_audio = audioread('dist_marshall_di.wav');
%     processed_audio = conv(input_audio, ir_band_filter);
%     processed_audio = processed_audio / max(abs(processed_audio));
%     audiowrite('processed_audio.wav', processed_audio, 48000);
%     return
%     
    

    subplot(2,1,1)
    plot(ir_band_filter)
    grid
    subplot(2,1,2)
    semilogx(20*log10(abs(fft(ir_band_filter, 48000))));
    grid


    return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    % 2 - Find  the non-linear characteristics: model of waveform compressor

    amount_of_periods = 500;

    % Find out sine frequency
    [~, sine_freq] = max(abs(fft(sine_input)));
    sine_freq = sine_freq - 1;

    sine_input = sine_input - mean(sine_input);
    sine_input = sine_input / max(abs(sine_input));
    sine_output = sine_output - mean(sine_output);
    sine_output = sine_output / max(abs(sine_output));


%     subplot(2,1,1)
%     plot(sine_input)
%     hold on
%     plot(sine_output)
%     grid

    % Align them on time
    phase_trim_init = 1;

    % Find the best range for the phase envelope, will be used for aligning
    % the signals on time by modifying the phase of the output
    best_variance = inf;
    best_trim_point = 2;
    for phase_trim_end = 2 : 1 : 20
        phase_fft_sine_in = unwrap(angle(fft(sine_input)));
        phase_fft_sine_in_trim = phase_fft_sine_in( round(length(phase_fft_sine_in) * phase_trim_init / 100) + 1 : round(length(phase_fft_sine_in) * phase_trim_end / 100));

        phase_fft_sine_out = unwrap(angle(fft(sine_output)));
        phase_fft_sine_out_trim = phase_fft_sine_out( round(length(phase_fft_sine_out) * phase_trim_init / 100) + 1 : round(length(phase_fft_sine_out) * phase_trim_end / 100));

        phase_envelope_in = [1 : length(phase_fft_sine_in_trim)]' \ phase_fft_sine_in_trim;
        phase_envelope_out = [1 : length(phase_fft_sine_out_trim)]' \ phase_fft_sine_out_trim;

        phase_envelope_difference = phase_envelope_out / phase_envelope_in;
        phase_fft_sine_out = phase_fft_sine_out / phase_envelope_difference;

        sine_output = real(ifft(abs(fft(sine_output)) .* [cos(phase_fft_sine_out) + 1i*sin(phase_fft_sine_out)]));
        sine_output = sine_output - mean(sine_output);

        % Check if the current case provides a better match between signals
        % on time
        new_correlation = (var(sine_output - sine_input));
        if new_correlation < best_variance;
            best_variance = new_correlation;
            best_trim_point = phase_trim_end;
        end
    end

    % Compute the phase of the spectrom of both signals
    phase_fft_sine_in = unwrap(angle(fft(sine_input)));
    phase_fft_sine_in_trim = phase_fft_sine_in( round(length(phase_fft_sine_in) * phase_trim_init / 100) + 1 : round(length(phase_fft_sine_in) * best_trim_point / 100));
    phase_fft_sine_out = unwrap(angle(fft(sine_output)));
    phase_fft_sine_out_trim = phase_fft_sine_out( round(length(phase_fft_sine_out) * phase_trim_init / 100) + 1 : round(length(phase_fft_sine_out) * best_trim_point / 100));

    % Find the envelope of both signals doing a linear regression
    phase_envelope_in = [1 : length(phase_fft_sine_in_trim)]' \ phase_fft_sine_in_trim;
    phase_envelope_out = [1 : length(phase_fft_sine_out_trim)]' \ phase_fft_sine_out_trim;

    % Correct the phase so they match
    phase_envelope_difference = phase_envelope_out / phase_envelope_in;
    phase_fft_sine_out = phase_fft_sine_out / phase_envelope_difference;

    % Recover output signal with corrected phase
    sine_output = real(ifft(abs(fft(sine_output)) .* [cos(phase_fft_sine_out) + 1i*sin(phase_fft_sine_out)]));
    sine_output = sine_output - mean(sine_output);

    % Compensate for a possible phase inversion in the alignment
    if corr(sine_input, sine_output) < 0
        sine_output = -sine_output;
    end

    % Trim to a part of them only
    sine_input = sine_input(round(length(sine_input)/4) : round(length(sine_input)/4) + round(sample_rate / sine_freq * amount_of_periods));
    sine_output = sine_output(round(length(sine_input)/4) : round(length(sine_input)/4) + round(sample_rate / sine_freq * amount_of_periods));

    % Eliminate the possible DC component
    sine_input = sine_input - mean(sine_input);
    sine_output = sine_output - mean(sine_output);
    % Normalize amplitudes
    sine_input = sine_input / max(abs(sine_input));
    sine_output = sine_output / max(abs(sine_output));


    
    
    
    
    % They need to be further vertically aligned
%     pp_input = max(sine_input) - min(sine_input);
%     pp_output = max(sine_output) - min(sine_output);
%     sine_input = sine_input /
    
%     central_value_input = (max(sine_input) + min(sine_input)) / 2;
%     central_value_output = (max(sine_output) + min(sine_output)) / 2;
%     sine_input = sine_input - central_value_input;
%     sine_output = sine_output - central_value_output;
    


%     close all
%     plot(sine_input(1000:1300));
%     hold on
%     plot(sine_output(1000:1300));
%     grid
%     return
    
    
    
    
    
    
    
    
    
    % Learn the transfer curve shape
    resolution = 2^16;

    % Scale both to the discrete range of values given by the resolution of
    % ADC and DAC
%     sine_input = round(sine_input * round(resolution/2));
%     sine_output = round(sine_output * round(resolution/2));

    % Average all periods
    points_to_match = 10;
    initial_points = sine_input(1 : points_to_match);
    current_transfer_curve = [];
    averaged_transfer_curve = [];
    count_averaged_curves = 0;
    counter = 1;
    for i = points_to_match : length(sine_input)

        if norm(sine_input(i-points_to_match+1 : i) - initial_points) < 0.1
            if ~isempty(current_transfer_curve)
                current_transfer_curve(counter, :) = [sine_input(i) sine_output(i)];
                counter = counter + 1;
                if ~isempty(averaged_transfer_curve)
                    if length(current_transfer_curve) == length(averaged_transfer_curve)
                        averaged_transfer_curve = averaged_transfer_curve + current_transfer_curve;
                        count_averaged_curves = count_averaged_curves + 1;
                    end
                else
                    averaged_transfer_curve = current_transfer_curve;
                    count_averaged_curves = 1;
                end
                current_transfer_curve = [];
            end
            counter = 1;
        else
            current_transfer_curve(counter, :) = [sine_input(i) sine_output(i)];
            counter = counter + 1;
        end
        
    end
    averaged_transfer_curve = averaged_transfer_curve / count_averaged_curves;



%     close all
%     plot([averaged_transfer_curve(:,1); averaged_transfer_curve(1,1)], [averaged_transfer_curve(:,2); averaged_transfer_curve(1,2)],'-', 'LineWidth', 1.5, 'Color', 'blue')
%     hold on
%     grid on
% return



    % Separate loop in two parts (up and down)
    low_point_radius = 0;
    high_point_radius = 0;
    for i = 1 : length(averaged_transfer_curve)
        current_radius = sqrt(averaged_transfer_curve(i, 1)^2 + averaged_transfer_curve(i, 2)^2);
        if current_radius > high_point_radius
            
            if averaged_transfer_curve(i, 1) > 0 && averaged_transfer_curve(i, 2) > 0
                high_point_radius = current_radius;
                hi_point_pos = i;
            end
        end
        if current_radius > low_point_radius
            if averaged_transfer_curve(i, 1) < 0 && averaged_transfer_curve(i, 2) < 0
                low_point_radius = current_radius;
                lo_point_pos = i;
            end
        end
    end
%     plot(averaged_transfer_curve(hi_point_pos, 1), averaged_transfer_curve(hi_point_pos, 2), 'x')
%     plot(averaged_transfer_curve(lo_point_hi, 1), averaged_transfer_curve(lo_point_hi, 2), '^')

    curve_part_1 = [];
    counter = 1;
    position = lo_point_pos;
    while position ~= hi_point_pos
        if position > length(averaged_transfer_curve)
            position = 1;
        end
        curve_part_1(counter, :) = averaged_transfer_curve(position, :);
        counter = counter + 1;
        position = position + 1;
    end
    curve_part_2 = [];
    counter = 1;
    position = hi_point_pos;
    while position ~= lo_point_pos
        if position > length(averaged_transfer_curve)
            position = 1;
        end
        curve_part_2(counter, :) = averaged_transfer_curve(position, :);
        counter = counter + 1;
        position = position + 1;
    end

    % Organize parts as rising and falling
    if curve_part_1(1) < curve_part_1(length(curve_part_1))
        rise_curve = curve_part_1;
        fall_curve = curve_part_2;
    else
        rise_curve = curve_part_2;
        fall_curve = curve_part_1;
    end


    % Find a rough intersection point
    cross_point_approx = [];
    for i = 1 : length(rise_curve)-1
        for j = 1 : length(fall_curve)-1
            [collision, collision_point, ~, ~] = detect_collision(rise_curve(i, :), rise_curve(i+1, :), ...
                                                 [fall_curve(j,:); fall_curve(j+1,:)], 0.001);
            if collision
                cross_point_approx = collision_point;
                cross_point_rise = i;
                cross_point_fall = j;
                break
            end
        end
    end

    % Find a fine intersection point
    [~, fine_collision_point, ~, ~] = detect_collision(rise_curve(cross_point_rise, :), rise_curve(cross_point_rise+1, :), ...
                                                 [fall_curve(cross_point_fall,:); fall_curve(cross_point_fall+1,:)], 0.00001);

    % Offset curves to make them cross on (0,0)
    for i = 1 : length(rise_curve)
        rise_curve(i,:)  = rise_curve(i,:) - fine_collision_point;
    end

    for i = 1 : length(fall_curve)
        fall_curve(i,:)  = fall_curve(i,:) - fine_collision_point;
    end

    % Normalize the amplitude to unity
    rise_curve = rise_curve / max(max(abs(rise_curve)));
    fall_curve = fall_curve / max(max(abs(fall_curve)));

    close all
    plot(rise_curve(:,1), rise_curve(:,2), 'r');
    hold on
    plot(fall_curve(:,1), fall_curve(:,2), 'b');

    grid on

%     [rise_curve_mix, fall_curve_mix] = scale_transfer_curves(rise_curve, fall_curve, 20);
%     plot(rise_curve_mix(:,1), rise_curve_mix(:,2), 'r', 'LineWidth', 1);
%     plot(fall_curve_mix(:,1), fall_curve_mix(:,2), 'b', 'LineWidth', 1);

    [rise_curve_mix, fall_curve_mix] = scale_transfer_curves(rise_curve, fall_curve, 50);
    plot(rise_curve_mix(:,1), rise_curve_mix(:,2), 'r', 'LineWidth', 1);
    plot(fall_curve_mix(:,1), fall_curve_mix(:,2), 'b', 'LineWidth', 1);

    [rise_curve_mix, fall_curve_mix] = scale_transfer_curves(rise_curve, fall_curve, 80);
    plot(rise_curve_mix(:,1), rise_curve_mix(:,2), 'r', 'LineWidth', 1);
    plot(fall_curve_mix(:,1), fall_curve_mix(:,2), 'b', 'LineWidth', 1);
    
    

    return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    % Find the central crossing
    minimum_difference = inf;
    for i = 1 : length(averaged_transfer_curve)
        for j = 1 : length(averaged_transfer_curve)
            if i ~= j
                if norm(averaged_transfer_curve(i) - averaged_transfer_curve(j)) < minimum_difference
                    minimum_difference = norm(averaged_transfer_curve(i) - averaged_transfer_curve(j));
                    central_point_1 = i;
                    central_point_2 = j;
                end
            end
        end
    end
    plot(averaged_transfer_curve(central_point_1, 1), averaged_transfer_curve(central_point_1, 2), 'x')
    plot(averaged_transfer_curve(central_point_2, 1), averaged_transfer_curve(central_point_2, 2), '^')


return






















    initial_point = sine_input(1);
    last_point_used = 1;
    while 1
        new_curve = [];
        previous_difference = inf;
        counter = 1;
        for i = last_point_used : length(sine_input)
            if (sine_input(i) - initial_point) < previous_difference
                previous_difference = sine_input(i) - initial_point;

                % Save new position
                 new_curve(counter, :) = [sine_input(i) sine_output(i)];
                 counter = counter + 1;

            else
                last_point_used = i;
                break;
            end
        end

        % Print new loop
        close all
        plot(new_curve(:,1), new_curve(:,2), 'x-');
        pause
    end

    

























    % Find start point
    for i = 1 : 1000
        [~, start_point] = min(abs(sine_input(i : 1000)));
        if sine_input(start_point + 1) > sine_input(start_point)
            central_pos_margin = 2 * sine_input(start_point);
            break
        end
    end

    close all
    plot(sine_input(1:1000));
    hold on
    plot(start_point, sine_input(start_point), 'x')
    grid
    return

















    % Find initial point and trim signals
    initial_point = 1;
    for i = 2 : length(sine_input)
        current_angle = atan2(sine_input(i), sine_input(i-1));
        if ( current_angle > 0 && current_angle < (pi/2)) && (sine_input(i) < 0.01)
            initial_point = i;
            break
        end
    end
    initial_point
    sine_input = sine_input(initial_point : length(sine_input));
    sine_output = sine_output(initial_point : length(sine_output));



    transfer_curve = zeros(length(sine_input), 2);

    for i = 1 : length(sine_input)

        x = sine_input(i);
        y = sine_output(i);
        transfer_curve(i,:) = [x, y];
%         plot(x,y, '.b')
%         hold on
%         transfer_curve(sine_input(i) + round(resolution/2) + 1) = transfer_curve(sine_input(i) + round(resolution/2) + 1) ...
%                                                                   + sine_output(i) + round(resolution/2) + 1;

    end

    close all
    plot(transfer_curve(:,1), transfer_curve(:,2))
    xlabel('Input value')
    ylabel('Output value')
    title('I/O transfer curve (distorted tube)')
    grid

%     hold on
%     plot([-1 1], [-1 1], 'LineWidth', 2);
    return
    
    
    
    
    
    
    
    
    
    
    
    
    


    subplot(2,1,2)
    plot(sine_input(4000:4800))
    hold on
    plot(sine_output(4000:4800))
    grid
    return

    
    plot(phase_fft_sine_in)
    hold on
    plot(phase_fft_sine_out)
    grid
    return
    
    

    plot(sine_input);
    grid
    return
    
    
    
    
    
    
    
    
    

end

function [collision, collision_point, reflection_direction, trazed_ray] = detect_collision(ray_origin, ray_end, boundary_limits, error_margin)

    % Find the boundary intermediate points
    boundary_direction = atan2(boundary_limits(2,2) - boundary_limits(1,2), boundary_limits(2,1) - boundary_limits(1,1));
    boundary_length = sqrt((boundary_limits(2,2) - boundary_limits(1,2))^2 + (boundary_limits(2,1) - boundary_limits(1,1))^2);
    boundary_points_step_distance = error_margin / 2;
    boundary_points = [boundary_limits(1,1) + cos(boundary_direction) * [0 : boundary_points_step_distance : boundary_length]' ...
                       boundary_limits(1,2) + sin(boundary_direction) * [0 : boundary_points_step_distance : boundary_length]' ];

    % Find the ray intermediate points
    ray_points_step_distance = error_margin / 2;
%     trazed_ray = [[ray_origin(1) : ray_points_step_distance : ray_end(1)]' ...
%                   [ray_origin(2) : ray_points_step_distance : ray_end(2)]'];

    ray_direction = atan2(ray_end(2) - ray_origin(2), ray_end(1) - ray_origin(1));
    max_radius = norm(ray_end(2) - ray_origin(2), ray_end(1) - ray_origin(1));
    trazed_ray = [ray_origin(1) + cos(ray_direction) * [0 : ray_points_step_distance : max_radius]' ...
                  ray_origin(2) + sin(ray_direction) * [0 : ray_points_step_distance : max_radius]'];

%     if sum(size(trazed_ray) == [1,1])
%         disp('');
%     end
%     size(boundary_points)
%     disp('---')
%               
    % Check if the ray gets close enough to any of the intermediate points of the boundary
    collision = 0;
    collision_point = [];
    reflection_direction = [];
    for i = [1 : size(trazed_ray, 1)]
        if min(sqrt( (trazed_ray(i, 1) - boundary_points(:,1)).^2 + (trazed_ray(i, 2) - boundary_points(:,2)).^2 )) <= error_margin
            collision = 1;

            % Old code: the collision point is approximate
            % collision_point = trazed_ray(i, :);

            % New code: the collision point is more exact
            [~, collision_point_index] = min(sqrt( (trazed_ray(i, 1) - boundary_points(:,1)).^2 + (trazed_ray(i, 2) - boundary_points(:,2)).^2 ));
            collision_point = boundary_points(collision_point_index, :);

            reflection_direction = 2 * boundary_direction - ray_end;
            return
        end        
    end

end


function [rise_mix, fall_mix] = scale_transfer_curves(rise_curve, fall_curve, scale_value)

    rise_mix = zeros(size(rise_curve));
    for i = 1 : length(rise_mix)
        rise_mix(i, :) = rise_curve(i,:) * scale_value/100 + [rise_curve(i,1) rise_curve(i,1)] * (100-scale_value) / 100;
    end
    rise_mix = rise_mix * scale_value / 100;

    fall_mix = zeros(size(fall_curve));
    for i = 1 : length(fall_mix)
        fall_mix(i, :) = fall_curve(i,:) * scale_value/100 + [fall_curve(i,1) fall_curve(i,1)] * (100-scale_value) / 100;
    end
    fall_mix = fall_mix * scale_value / 100;

end
