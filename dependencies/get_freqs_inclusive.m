function freqs_incl = get_freqs_inclusive(freqs, freq_rng)
    freqs_incl = get_freqs_exclusive(freqs, freq_rng);
    if mod(freq_rng(1), freqs(2))
        freqs_incl(find(freqs < freq_rng(1), 1, 'last')) = 1;
    end
    if mod(freq_rng(2), freqs(2))
        freqs_incl(find(freqs > freq_rng(2), 1, 'first')) = 1;
    end
end
