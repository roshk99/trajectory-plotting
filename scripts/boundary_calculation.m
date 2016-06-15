function values = boundary_calculation(smooth_data, set_num, ...
    canalview, howmany, idx1, idx2, plotSurface)
    
    threshold = 0.01;
    plotBoundaries = false;
    plotSurface_internal = false;
    filterSurface_internal = false;
    
    numDemos = numel(howmany);
    nPoints = size(smooth_data{1}, 1);
    
    allXs = zeros(nPoints,numDemos);
    allYs = zeros(nPoints,numDemos);
    allZs = zeros(nPoints,numDemos);
    for ii = 1:numDemos
        allXs(:,ii) = smooth_data{ii}(:,1);
        allYs(:,ii) = smooth_data{ii}(:,2);
        allZs(:,ii) = smooth_data{ii}(:,3);
    end
    
    B = calculate_mean(allXs, allYs, allZs);
    
    [Router, xyz_distance] = find_boundaries([B(2).xmean, ...
        B(2).ymean, B(2).zmean], allXs, allYs, allZs);
    
    if plotBoundaries
        plot_boundaries(data, B, xyz_distance, numDemos, threshold);
    end
    
    %Calculate and plot canal surface
    t = linspace(0,1,nPoints).';
    tt = t(idx1:end-idx2);
    xx2 = B(2).xmean(idx1:end-idx2);
    yy2 = B(2).ymean(idx1:end-idx2);
    zz2 = B(2).zmean(idx1:end-idx2);
    RRouter = Router(idx1:end-idx2);
    
    [~,canalouter,T,~,~,N2,B2] = ...
        canalSurface(tt,xx2,yy2,zz2,RRouter,plotSurface_internal,...
        filterSurface_internal);
    
    if plotSurface
        plot_surface(xx2, yy2, zz2, smooth_data, canalouter, ...
            canalview, set_num, numDemos)
    end
    
    values = struct([]);
    values(1).xmean = xx2;
    values(1).ymean = yy2;
    values(1).zmean = zz2;
    values(1).Router = Router(idx1:end-idx2);
    values(1).xyz_distance = xyz_distance;
    values(1).N2 = N2';
    values(1).B2 = B2';
    values(1).T = T';
end

function B = calculate_mean(allXs, allYs, allZs)
    B = struct([]);
    B(2).xmean = mean(allXs,2);
    B(2).ymean = mean(allYs,2);
    B(2).zmean = mean(allZs,2);
end

function plot_boundaries(data, B, xyz_distance, numDemos, threshold)
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

function plot_surface(xx2, yy2, zz2, data, canalouter, canalview, ...
    set_num, numDemos)
    
    figure;hold on
    for kk = 1:size(canalouter,3)
        C = canalouter(:,:,kk);
        plot3(C(1,:),C(2,:),C(3,:),'k');
    end
    plot3(xx2,yy2,zz2,'b','linewidth',2);
%     surf(squeeze(canalouter(1,:,:)), squeeze(canalouter(2,:,:)), squeeze(canalouter(3,:,:)));
%     shading interp;
%     alpha 0.5;
%     colormap bone;
    
    for ii=1:numDemos
        plot3(data{ii}(:,1), data{ii}(:,2), data{ii}(:,3),'r',...
            'linewidth',2);
    end
    title(sprintf('Set %i - Canal Surface', set_num));
    axis equal
    view(canalview);
    set(gcf, 'Position', get(0, 'Screensize'));
    hold off;
end
