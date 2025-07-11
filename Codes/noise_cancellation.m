function varargout = noise_cancellation(noisy_speech, external_noise, clean_speech, tonal_freqs, mode)
    fs = 44.1e3;
    n1 = length(external_noise); % length of the vector
    lambda = 0.9995; % forgetting factor
    p = 11; % order
    N = p+1; % number of taps
    delta = 10^-4;

    % initialization
    w_n = zeros(N, 1);

    % Initializing P_n to 1 / delta * I
    P_n = 1 / delta * eye(N); 

    % g_n = zeros(N, 1);
    v_hat = zeros(n1, 1);
    filtered_speech = zeros(n1, 1);

    if lower(mode) == "partial"
        n_freqs = length(tonal_freqs);
        frequencies = tonal_freqs / fs;
        w = 2 * pi * frequencies;
        r = 0.999;
        b = ones(n_freqs, 3);
        a = ones(n_freqs, 3);
        for j = 1:n_freqs
            b(j, 2) = -2 * cos(w(j));
            a(j, 2) = -2 * r * cos(w(j));
            a(j, 3) = r^2;
        end
        outputs_notched = zeros(n_freqs, n1);
    end

    for i = 1 : n1
        x_ = zeros(N, 1);
        if lower(mode) == "partial" 
            if i == 1
                o_1 = 0;
                o_2 = 0;
                i_0 = external_noise(1);
                i_1 = 0;
                i_2 = 0;
            elseif i == 2
                o_1 = outputs_notched(1, 1);
                o_2 = 0;
                i_0 = external_noise(2);
                i_1 = external_noise(1);
                i_2 = 0;
            else
                o_1 = outputs_notched(1, i-1);
                o_2 = outputs_notched(1, i-2);
                i_0 = external_noise(i);
                i_1 = external_noise(i-1);
                i_2 = external_noise(i-2);
            end
            outputs_notched(1, i) = 1/a(1, 1) * (-a(1, 2:3) * [o_1; o_2] + b(1, 1:3) * [i_0; i_1; i_2]);
            for k = 2:n_freqs
                if i == 1
                    o_1_k = 0;
                    o_2_k = 0;
                    i_0_k = outputs_notched(k-1, 1);
                    i_1_k = 0;
                    i_2_k = 0;
                elseif i == 2
                    o_1_k = outputs_notched(k, 1);
                    o_2_k = 0;
                    i_0_k = outputs_notched(k-1, 2);
                    i_1_k = outputs_notched(k-1, 1);
                    i_2_k = 0;
                else
                    o_1_k = outputs_notched(k, i-1);
                    o_2_k = outputs_notched(k, i-2);
                    i_0_k = outputs_notched(k-1, i);
                    i_1_k = outputs_notched(k-1, i-1);
                    i_2_k = outputs_notched(k-1, i-2);
                end
                outputs_notched(k, i) = 1/a(k, 1) * (- a(k, 2:3) * [o_1_k; o_2_k] + b(k, 1:3) * [i_0_k; i_1_k; i_2_k]); 
            end
            
            k = 1;
            for j = i : -1 : i - N + 1
                if j <= 0
                    x_(k) = 0;
                    k = k + 1;
                else
                    x_(k) = outputs_notched(n_freqs, j);
                    k = k + 1;
                end
            end
        elseif lower(mode) == "full"
            k = 1;
            for j = i : -1 : i - N + 1
                if j <= 0
                    x_(k) = 0;
                    k = k + 1;
                else
                    x_(k) = external_noise(j);
                    k = k + 1;
                end
            end
        end

        v_hat(i) = x_' * w_n;

        alpha = noisy_speech(i) - v_hat(i);
        filtered_speech(i) = alpha;

        g_n = (P_n * x_) * 1 / (lambda + x_' * P_n * x_);

        P_n = (1 / lambda) * ( P_n - g_n * x_' * P_n);

        w_n = w_n + alpha * g_n;
    end
    if nargout >= 1
        varargout{1} = filtered_speech;
    end
    if nargout == 4
        % calculates the notched counterparts of the signals
        clean_speech_notched = clean_speech;
        noisy_speech_notched = noisy_speech;
        filtered_speech_notched = filtered_speech;
        for i = 1:n_freqs
            clean_speech_notched = filter(b(i, 1:3), a(i, 1:3), clean_speech_notched);
            noisy_speech_notched = filter(b(i, 1:3), a(i, 1:3), noisy_speech_notched);
            filtered_speech_notched = filter(b(i, 1:3), a(i, 1:3), filtered_speech_notched);
        end 
        varargout{2} = clean_speech_notched;
        varargout{3} = noisy_speech_notched;
        varargout{4} = filtered_speech_notched;
    end
end