function Button_Down_Fcn_Bar(BarHandles, Event)
    % Determine which bar was clicked
    idx = round(Event.IntersectionPoint(1));

    % Adjust bar color
    BarHandles.CData(:, :) = repmat([0, 1, 0], [size(BarHandles.CData, 1), 1]);
    BarHandles.CData(idx, :) = [1, 0, 0];
    BarHandles.FaceColor = 'flat';

    % Adjust label color
    COLOR_STR = '\color{red}';
    BarHandles.Parent.XTickLabel = cellfun(@(a) strrep(a, COLOR_STR, ''), ...
        BarHandles.Parent.XTickLabel, 'UniformOutput', false);
    copy_str = BarHandles.Parent.XTickLabel{idx};
    BarHandles.Parent.XTickLabel{idx} = ...
        [COLOR_STR, BarHandles.Parent.XTickLabel{idx}];

    % Store respective x-tick-label and y-value in clipboard
    copy_str = strrep(copy_str, '\_', '_');
    copy_str = sprintf('%s %f', copy_str, BarHandles.YData(idx));
    clipboard('copy', copy_str);
end
