function freqs_excl = get_freqs_exclusive(freqs, freq_rng)
    if isempty(freq_rng)
        freqs_excl = [];
        return;
    end
    assert(length(freq_rng) == 2);

    freqs_excl = freq_rng(1) <= freqs & freqs <= freq_rng(2);
end
