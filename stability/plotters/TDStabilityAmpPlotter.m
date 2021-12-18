function [figure_handle] = TDStabilityRegionPlotter(data_raw, z1_r, z2_r, options)
%TDSTABILITYPLOTTER produces the plot for a two-dimensional stability plot
%   data_raw (struct) - data calculated using OverlayedStabilityData
%   z1_r  (vector) - radius of z1 component. 
%   z2_r  (vector) - radius of z2 component. 

if(nargin == 3)
    options = struct();
end
default_options = {
    {'FigureIndex',             1}
    {'ClearFigure',             true}
    {'DrawAxis',                true}
    {'Overlay',                 'all'}
    {'OverlayEdgeColor',        .3 * [1 1 1]}
    {'OverlayLineWidth',        1}
    {'Union',                   'all'}
    {'UnionEdgeColor',          .3 * [1 1 1]}
    {'UnionBackgroundColor',    .8 * [1 1 1]} 
    {'UnionLineWidth',          1}
    {'IndexedOverlayEdgeStyle', {}}  % nx2 cell array A{j,:} = {z1_special, line_style}. contour with z1=z1_special will be created using arguments line_style{:}
    {'ContourLevels'           [0 1]}
    {'FontSize',               12}
    {'FontName',               'Minion Pro'}
    {'AxisSquare',             true}
    {'NumTicks',               5}
};
options = setDefaultOptions(options, default_options);


% -- Set Figure Index --------------------------------------------------------------------------------------------------
if(isempty(options.FigureIndex))
    figure_handle = figure();
else
    figure_handle = figure(options.FigureIndex);
    if(options.ClearFigure)
        clf;
    end    
end
hold on;

% -- draw axis ---------------------------------------------------------------------------------------------------------
if(options.DrawAxis)
    plot([min(z1_r) max(z1_r)], [0 0], 'k'); 
    plot([0 0], [min(z2_r) max(z2_r)], 'k');
end

% -- create overlay plot -----------------------------------------------------------------------------------------------
if(min(data_raw(:)) < max(options.ContourLevels)) % only plot if at least one value is below threshold
    
%     surf(z1_r, z2_r, log(data_raw)); shading interp; hold on;
%     contour(z1_r, z2_r, data_raw, options.ContourLevels, 'color', [0 0 0]); hold off;
%     colorbar; caxis([-5,5]);
%     view([30 25]);
%     zlim([-5,5]);
    
    
    contourf(z1_r, z2_r, data_raw, options.ContourLevels, 'color', options.UnionEdgeColor, 'linewidth', options.UnionLineWidth);
    if(~isempty(options.UnionBackgroundColor))
        colormap([options.UnionBackgroundColor; 1 1 1]);
    end
else
    warning('there were not stable points in union set');
end

% -- set Font and axis ----
set(gca, 'FontSize', options.FontSize, 'FontName', options.FontName);

if(options.AxisSquare)
    axis square;
end

xtick_vec = linspace(min(z1_r), max(z1_r), options.NumTicks);
ytick_vec = linspace(min(z2_r), max(z2_r), options.NumTicks);
xticks(xtick_vec);
yticks(ytick_vec);
xlabel('$|k_1|$', 'interpreter', 'latex');
ylabel('$|k_2|$', 'interpreter', 'latex');

hold off;
end