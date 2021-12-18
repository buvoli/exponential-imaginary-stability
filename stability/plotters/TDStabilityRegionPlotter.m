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
    {'ThreeDimensionalAmp'     false}
    {'ZLim'                    []}
    {'ShowUnstable'            false}
    {'ShowColorbar',           true}
    {'LabelAxis',              true}
    {'AmpCaxis',               [-5,0]}
    {'LogAmp',                 true}
    {'View',                   [0 90]}
    {'ColorMap',               []}
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
    if(options.ThreeDimensionalAmp)
        amp_to_plot = data_raw;
        if( ~ options.ShowUnstable)
            amp_to_plot(amp_to_plot > 1) = NaN;
        end
        if( options.LogAmp )
            amp_to_plot = log(amp_to_plot);
            clz = 0;
        else
            clz = 1;
        end
        surf(z1_r, z2_r, amp_to_plot); shading interp; hold on;
        [~, cnt] = contour(z1_r, z2_r, data_raw, options.ContourLevels, 'color', [0 0 0]); hold off;
        if(isprop(cnt, 'ZLocation'))
            cnt.ZLocation = clz;
        elseif(isprop(cnt, 'ContourZLevel'))
            cnt.ContourZLevel = clz;
        end
        view(options.View);
        if(options.ShowColorbar)
            colorbar; 
        end
        if(options.ZLim)
            zlim(options.ZLim);
        end
        caxis(options.AmpCaxis);
    else
        data_raw_trim = data_raw;
        data_raw_trim(data_raw > max(options.ContourLevels)) = NaN;
        contourf(z1_r, z2_r, data_raw_trim, options.ContourLevels, 'color', options.UnionEdgeColor, 'linewidth', options.UnionLineWidth); hold on;
        contour(z1_r, z2_r, data_raw, options.ContourLevels, 'color', options.UnionEdgeColor, 'linewidth', options.UnionLineWidth); hold off;
        if(~isempty(options.ColorMap))
            colormap([options.ColorMap])
        elseif(~isempty(options.UnionBackgroundColor))
            colormap([options.UnionBackgroundColor; 1 1 1]);
        end
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

if(options.LabelAxis)
    xlabel('$k_1$', 'interpreter', 'latex');
    ylabel('$k_2$', 'interpreter', 'latex');
    if(options.ThreeDimensionalAmp)
        zlabel('$|R(ik_1,ik_2)|$', 'interpreter', 'latex');
    end
end

hold off;
end