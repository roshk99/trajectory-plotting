function plotting(raw_data, data, mean_vals, canal, traj, Router, ...
    plotRawTrajectories, plotSmoothTrajectories, plotComponents, ...
    plotSurface, plotNewTrajectories, plotBoundaries, how_many, ...
    plot_method, handles)

new_data = {};
new_data_raw = {};
counter = 1;
for ind1=how_many
    new_data{counter} = data{ind1};
    new_data_raw{counter} = raw_data{ind1};
    counter = counter+1;
end
data = new_data;
raw_data = new_data_raw;

if plotRawTrajectories
    plot_data(raw_data, view_vec, plotComponents);
end

if plotSmoothTrajectories
    plot_data(data, view_vec, plotComponents);
end

if plotSurface
    plot_surface(mean_vals, data, canal, plot_method, handles);
end

if plotNewTrajectories
    plot_trajectories(traj);
end

if plotBoundaries
    plot_boundaries(data, mean_vals, Router, traj);
end
end

function plot_data(data, plotComponents)
%Plot the demonstrations and the demonstrations separated into components
cc = jet(length(data));

figure; hold on;
legend_vec = {};
for ind1 = 1:length(data)
    plot3(data{ind1}(:, 1), data{ind1}(:,2), data{ind1}(:,3), ...
        'color', cc(ind1,:), 'linewidth', 1.5);
    legend_vec{ind1} = sprintf('Series %i', ind1);
end
hold off;
legend(legend_vec);
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
axis equal;
set(gcf, 'Position', get(0, 'Screensize'));

if plotComponents
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(3,1,1); hold all;
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

function plot_boundaries(data, mean_vals, Router, traj)
%Plots the boundaries in x, y, and z direction (only for circles)

figure;
subplot(3,1,1);
hold on;
for ii=1:length(data)
    plot(data{ii}(:,1),'r');
end
plot(mean_vals{1},'b','Linewidth', 1.5);
for ii=1:length(traj)
   plot(traj{ii}(:,1),'g'); 
end
threshold = 0.01;
plot(mean_vals{1}+Router+threshold, 'k', 'Linewidth', 1.5);
plot(mean_vals{1}-Router-threshold, 'k', 'Linewidth', 1.5);
ylabel('X');axis tight; grid on; hold off;

subplot(3,1,2);
hold on;
for ii=1:length(data)
    plot(data{ii}(:,2),'r');
end
plot(mean_vals{2},'b', 'Linewidth', 1.5);
for ii=1:length(traj)
   plot(traj{ii}(:,2),'g'); 
end
threshold = 0.01;
plot(mean_vals{2}+Router+threshold, 'k', 'Linewidth', 1.5);
plot(mean_vals{2}-Router-threshold, 'k', 'Linewidth', 1.5);
ylabel('Y');axis tight;grid on; hold off;

subplot(3,1,3);
hold on;
for ii=1:length(data)
    plot(data{ii}(:,3),'r');
end
plot(mean_vals{3},'b', 'Linewidth', 1.5);
for ii=1:length(traj)
   plot(traj{ii}(:,3),'g'); 
end
threshold = 0.01;
plot(mean_vals{3}+Router+threshold, 'k', 'Linewidth', 1.5);
plot(mean_vals{3}-Router-threshold, 'k', 'Linewidth', 1.5);
ylabel('Z');axis tight;grid on; box on; hold off;
end

function plot_surface(mean_vals, data, canal, plot_method, handles)
%Plots the canal surface, directrix, and demonstrations

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
    set(gcf, 'Position', get(0, 'Screensize'));
end

hold off;
end

function plot_trajectories(trajs)
%Plots the reproductions
hold on;
for ii = 1:length(trajs)
    plot3(trajs{ii}(:,1), trajs{ii}(:,2), trajs{ii}(:,3), 'g', ...
        'Linewidth', 2);
end
hold off;
end