function strs = get_SH_strings(N_max)
    strs = {};
    for n = 0 : N_max
        for m = -n : n
            strs = [strs, sprintf('[%d,%+d]', n, m)]; %#ok<AGROW> 
        end
    end
end
