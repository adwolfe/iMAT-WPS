function h = plotReactionSet(height, ub_PFD, lb_PFD, ub_OFD, lb_OFD, mytext)
% Create the horizontal bar plot
h = figure;

num_plots = size(height,1); % Number of plots

% We need to adjust the bar and the errorbar calls to loop over each group
colors = [255, 111, 0;
          255,165,0;
          255,165,0;
          255, 227, 177;
          180 180 180;]/255; % Colors for the groups
for p = 1:num_plots % Loop over plots
    subplot(num_plots, 1, p);
    hold on;
    % xline(0,'-b', 'LineWidth',2)

    for i = 1:size(height, 2)
        fill( [lb_PFD(p, i); ub_PFD(p, i); ub_PFD(p, i); lb_PFD(p, i)], [i-0.4 i-0.4 i+0.4 i+0.4], colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
        fill( [lb_OFD(p, i); ub_OFD(p, i); ub_OFD(p, i); lb_OFD(p, i)], [i-0.4 i-0.4 i+0.4 i+0.4], nan, 'EdgeColor', [0,0,0],'LineWidth', 0.5);
        line([height(p, i) height(p, i)], [i-0.4 i+0.4],   'Color', 'r', 'LineWidth', 2); % Line from lower bound to height        
    end
    text(0, 5.5, '*', 'Color', 'blue', 'HorizontalAlignment', 'center','FontSize',25)

    hold off;
    title(mytext(p))

    xvals1 = [lb_PFD(p, 1:(end-1)), ub_PFD(p, 1:(end-1))]; % exclude the PFD of base bc it's almost unbounded
    xvals2 = [lb_OFD(p, 1:(end-1)), ub_OFD(p, 1:(end-1))]; % also for OFD
    mymin = min([0, max([min(xvals2) - abs(min(xvals2))*0.5, min(xvals1)])]);
    mymax = max([0, min([max(xvals2) + abs(max(xvals2))*0.5, max(xvals1)])]);

    xlim([mymin, mymax])
    set(gca,'YTickLabel',{'Exp+resp+sim','Exp+resp','Exp+sim','Exp','Biomass'},'YTick',[1 2 3 4 5])

end
xlabel('flux')
