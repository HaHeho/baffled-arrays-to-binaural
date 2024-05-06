function [ticks, ticklabels] = get_freqs_ticks(add_tick_freqs)

    if nargin < 1; add_tick_freqs = []; end

    ticks = unique([1, 4, 10, 40, 100, 400, 1e3, 4e3, 10e3, 20e3, 40e3, add_tick_freqs], 'sorted');
    if nargout > 1
        % transform into cell of strings
        ticklabels = arrayfun(@(tick) num2str(tick, '%.0f'), ticks, 'UniformOutput', false);

        % replace thousand abbreviations from the back
        ticklabels = cellfun(@(tick) flip(replace(flip(tick), '000', 'k')), ticklabels, 'UniformOutput', false);
    end

end
