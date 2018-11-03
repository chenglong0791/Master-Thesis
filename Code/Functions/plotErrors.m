function [] = plotErrors(t, errors, rmsErrors, variances, legendNames, titleName, kidnapIndex)
%PLOTERRORS Plots the error and the root mean-squared error

figure()
hold on;

for i = 1:size(errors, 1)
    
    if sum(errors(i, :)) == 0
        continue;
    end
    
    handle = plot(t, errors(i, :), 'LineWidth', 1.5, 'DisplayName', ['Error ', char(legendNames(i))]);
    set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex') - 1);
    
    if ~isempty(variances)
        fill([t, fliplr(t)],[errors(i, :) - variances(i, :), fliplr(errors(i, :) + variances(i, :))], ...
            handle.Color, 'FaceAlpha', 0.2, 'LineStyle','none', 'DisplayName', ['Variance ', char(legendNames(i))]);
        set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex') - 1);
    end
    
    %     plot([t(1), t(end)], [rmsErrors(i, :), rmsErrors(i, :)], '--', 'LineWidth', 1, ...
    %         'Color', handle.Color, 'DisplayName', ['RMSE ', char(legendNames(i))]);
    set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex') + 1);
end

hold off;

yl = ylim;

if ~isempty(kidnapIndex)
    line([kidnapIndex, kidnapIndex], [yl(1) yl(2)], 'Color', [0.5 0.5 0.5], 'LineStyle','--', 'DisplayName', 'Kidnapping')
end

title(titleName);
xlabel('Time in s');
ylabel('Error in m');
legend('-DynamicLegend');

figure()
hold on;

for i = 1:size(errors, 1)
    
    if sum(errors(i, :)) == 0
        continue;
    end
    
    handle = plot(t, errors(i, :), 'LineWidth', 1.5, 'DisplayName', ['Error ', char(legendNames(i))]);
    set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex') - 1);
    
    %     plot([t(1), t(end)], [rmsErrors(i, :), rmsErrors(i, :)], '--', 'LineWidth', 1, ...
    %         'Color', handle.Color, 'DisplayName', ['RMSE ', char(legendNames(i))]);
    set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex') + 1);
    
    
end

set(gca, 'YScale', 'log');
set(gca, 'YGrid', 'on')

hold off;

yl = ylim;

if ~isempty(kidnapIndex)
    line([kidnapIndex, kidnapIndex], [yl(1) yl(2)], 'Color', [0.5 0.5 0.5], 'LineStyle','--', 'DisplayName', 'Kidnapping')
end

title([titleName, ' (logarithmic scale)']);
xlabel('Time in s');
ylabel('Error in m');
legend('-DynamicLegend');

end


