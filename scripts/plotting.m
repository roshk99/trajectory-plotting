function plotting(raw_data, data, mean_vals, canal, traj, ...
    plotRawTrajectories, plotSmoothTrajectories, plotComponents, ...
    plotSurface, plotNewTrajectories, plot_method, set_num, ...
    fit_type, view_vec, handles)


if plotRawTrajectories
    plotTitle1 = sprintf('Series %i - Raw Data', set_num);
    plot_data(raw_data, view_vec, plotComponents, plotTitle1);
end

if plotSmoothTrajectories
    plotTitle1 = sprintf('Series %i - Smoothed Data', set_num);
    plot_data(data, view_vec, plotComponents, plotTitle1);
end


%Plot the canal surface
plotTitle2 = sprintf('Series %i using %s', set_num, fit_type);

if plotSurface
    plot_surface(mean_vals, data, canal, view_vec, ...
        plot_method, handles, plotTitle2);
end

if plotNewTrajectories
    plot_trajectories(traj);
end

end

function plot_data(data, view_vec, plotComponents, plotTitle)

cc = jet(length(data));

figure; hold on;
legend_vec = {};
for ind1 = 1:length(data)
    plot3(data{ind1}(:, 1), data{ind1}(:,2), data{ind1}(:,3), ...
        'color', cc(ind1,:), 'linewidth', 1.5);
    legend_vec{ind1} = sprintf('Series %i', ind1);
end
hold off;
title(sprintf('%s - 3D Trajectories', plotTitle));
legend(legend_vec);
xlabel('X'); ylabel('Y'); zlabel('Z');
view(view_vec);
grid on;
axis equal;
set(gcf, 'Position', get(0, 'Screensize'));

if plotComponents
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(3,1,1); hold all; title(sprintf('%s - X Component', plotTitle));
    subplot(3,1,2); hold all; title('Y Component');
    subplot(3,1,3); hold all; title('Z Component');
    
    for ind1 = 1:length(data)
        tvec = -(length(data{ind1}(:,1))-1):0;
        subplot(3,1,1);
        plot(tvec, data{ind1}(:,1), 'color', cc(ind1,:), ...
            'linewidth', 1.2);
        subplot(3,1,2);
        plot(tvec, data{ind1}(:,2), 'color', cc(ind1,:), ...
            'linewidth', 1.2);
        subplot(3,1,3);
        plot(tvec, data{ind1}(:,3), 'color', cc(ind1,:), ...
            'linewidth', 1.2);
    end
    legend(legend_vec);
    subplot(3,1,1); hold off; axis tight;
    subplot(3,1,2); hold off; axis tight;
    subplot(3,1,3); hold off; axis tight;
end

end

function plot_boundaries(data, B, xyz_distance, numDemos, threshold)
%Plots the boundaries in x, y, and z direction

%Plots the mean, trajectories, and boundary with an extra threshold
figure;
subplot(2,3,1);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,1),'k');
end
plot(B(2).xmean, 'b');
plot(xyz_distance(:,1)+B(2).xmean + threshold, 'r');
plot(-xyz_distance(:,1)+B(2).xmean - threshold, 'r');
title(sprintf('Set %i - X', set_num));axis square;grid;
hold off;

subplot(2,3,2);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,2),'k');
end
plot(B(2).ymean, 'b');
plot(xyz_distance(:,2)+B(2).ymean + threshold, 'r');
plot(-xyz_distance(:,2)+B(2).ymean - threshold, 'r');
title(sprintf('Set %i - Y', set_num));axis square;grid;
hold off;

subplot(2,3,3);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,3),'k');
end
plot(B(2).zmean, 'b');
plot(xyz_distance(:,3)+B(2).zmean + threshold, 'r');
plot(-xyz_distance(:,3)+B(2).zmean - threshold, 'r');
title(sprintf('Set %i - Z', set_num));axis square;grid;
hold off;

%Plots the data - mean and boundary
subplot(2,3,4);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,1)-B(2).xmean,'k');
end
plot(xyz_distance(:,1) + threshold, 'b');
title('Bound X');axis square;grid;
hold off;

subplot(2,3,5);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,2)-B(2).xmean,'k');
end
plot(xyz_distance(:,2) + threshold, 'b');
title('Bound Y');axis square;grid;
hold off;

subplot(2,3,6);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,3)-B(2).xmean,'k');
end
plot(xyz_distance(:,1) + threshold, 'b');
title('Bound Z');axis square;grid;
hold off;
end

function plot_surface(mean_vals, data, canal, view_vec, ...
    plot_method, handles, plotTitle)
%Plots the canal surface

if isa(handles, 'struct')
    axes(handles.axes1);
    cla reset;
else
    figure;
end

hold on;
if strcmp('circles', plot_method)
    for kk = 1:size(canal,3)
        C = canal(:,:,kk);
        plot3(C(1,:),C(2,:),C(3,:),'k');
    end
else
    surf(squeeze(canal(1,:,:)), squeeze(canal(2,:,:)), ...
        squeeze(canal(3,:,:)));
    shading interp;
    alpha 0.5;
    colormap bone;
end

plot3(mean_vals{1},mean_vals{2},mean_vals{3},'b','linewidth',2);

for ii=1:length(data)
    plot3(data{ii}(:,1), data{ii}(:,2), data{ii}(:,3),'r',...
        'linewidth',2);
end

if ~isa(handles, 'struct')
    title(plotTitle);
    set(gcf, 'Position', get(0, 'Screensize'));
end

xlabel('X'); ylabel('Y'); zlabel('Z');
view(view_vec);
axis equal;
hold off;
end

function plot_trajectories(trajs)
hold on;
for ii = 1:length(trajs)
    plot3(trajs{ii}(:,1), trajs{ii}(:,2), trajs{ii}(:,3), 'g', ...
        'Linewidth', 2);
end
hold off;
end