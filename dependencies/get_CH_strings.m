function strs = get_CH_strings(N_max)
    strs = {sprintf('[%d,%+d]', 0, 0)};
    for n = 1 : N_max
        strs = [strs, sprintf('[%d,%+d]', n, -n), sprintf('[%d,%+d]', n, n)]; %#ok<AGROW>
    end
end
