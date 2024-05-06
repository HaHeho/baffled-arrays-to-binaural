function irs = get_SOFA_IRs(sofa, idx)
    if nargin < 2; idx = 1 : size(sofa.Data.IR, 1); end

    irs = squeeze(permute(sofa.Data.IR(idx, :, :), [3, 1, 2]));
end
