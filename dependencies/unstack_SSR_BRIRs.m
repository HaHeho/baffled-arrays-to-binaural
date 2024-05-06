function sig = unstack_SSR_BRIRs(sig)
    sig = permute(reshape(sig, size(sig, 1), 2, []), [1, 3, 2]);
end
