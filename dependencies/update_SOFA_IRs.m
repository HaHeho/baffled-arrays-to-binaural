function sofa = update_SOFA_IRs(sofa, irs)
    sofa.Data.IR = permute(irs, [2, 3, 1]);
end
