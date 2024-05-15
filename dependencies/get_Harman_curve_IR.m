function ir = get_Harman_curve_IR(ir_len, fs, file)

    if nargin < 3
        file = 'resources/HRIR_equalizations/Harman_48k.wav';
    end

    is_even = ~mod(ir_len, 2);
    target_freqs = linspace(0, fs/2, ir_len/2 + 1).';

    if endsWith(file, 'txt')
        data = readmatrix(file);
        freqs = data(:, 1);
        spec_dB = data(:, 2);

    elseif endsWith(file, 'mat')
        data = load(file);
        spec_dB = AKboth2singleSidedSpectrum(db(abs(fft(data.IR))));
        freqs = linspace(0, data.srate/2, length(spec_dB));

    elseif endsWith(file, 'wav')
        [data, data_fs] = audioread(file);
        spec_dB = AKboth2singleSidedSpectrum(db(abs(fft(data))));
        freqs = linspace(0, data_fs/2, length(spec_dB));

    else
        error('Unknown file type "%s".', file);
    end

    % Interpolate one-sided magnitude spectrum to target frequency sampling
    target_spec_dB = interp1(freqs, spec_dB, target_freqs, 'nearest', 'extrap');
    % target_spec_dB = interp1(freqs, spec_dB, target_freqs, 'linear', 'extrap');
    target_spec_dB(1) = 0; % set 0 Hz bin
    % AKf;
    % semilogx(target_freqs, target_spec_dB);
    % hold on;
    % semilogx(freqs, spec_dB);
    % drawnow;

    ir = ifft(AKsingle2bothSidedSpectrum(db2mag(target_spec_dB), is_even));
    % AKf;
    % % AKp(data.IR, {'et2d'; 'm2d'}, 'fs', data.srate, 'c', 'r');
    % AKp(data, {'et2d'; 'm2d'}, 'fs', data_fs, 'c', 'r');
    % AKp(ir, {'et2d'; 'm2d'}, 'fs', fs);
    % legend({'Input', 'Interpolated'});
    % drawnow;
end
