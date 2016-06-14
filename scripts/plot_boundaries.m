function plot_boundaries(smooth_data, set_num, canalview)
    plotIt = false;
    filterIt = false;

    %% Calculate mean
    howmany = [1 2 3 4];
    
    numDemos = numel(howmany); 
    D = struct([]);
    for ii = 1:numDemos
        D(ii).Xs = smooth_data{ii}(:,1);
        D(ii).Ys = smooth_data{ii}(:,2);
        D(ii).Zs = smooth_data{ii}(:,3);
    end
    nPoints = numel(D(1).Xs);
    allXs = zeros(nPoints,numDemos);
    allYs = zeros(nPoints,numDemos);
    allZs = zeros(nPoints,numDemos);
    for ii=1:numDemos
        allXs(:,ii) = D(ii).Xs;
        allYs(:,ii) = D(ii).Ys;
        allZs(:,ii) = D(ii).Zs;
    end

    B = struct([]);
    B(2).xmean = mean(allXs,2);
    B(2).ymean = mean(allYs,2);
    B(2).zmean = mean(allZs,2);
%     rx = B(2).xmean.';
%     ry = B(2).ymean.';
%     rz = B(2).zmean.';

    %rr = [rx;ry;rz];
    
    threshold = 0.01;
    [Router, xyz_distance] = find_boundaries([B(2).xmean, ...
        B(2).ymean, B(2).zmean], allXs, allYs, allZs);
    
    %% plot the boundaries
    figure;
    subplot(2,3,1);
    hold on;
    for ii=1:numDemos
        plot(D(ii).Xs,'k');
    end
    plot(B(2).xmean, 'b');
    plot(xyz_distance(:,1)+B(2).xmean + threshold, 'r');
    plot(-xyz_distance(:,1)+B(2).xmean - threshold, 'r');
    title(sprintf('Set %i - X', set_num));axis square;grid;
    hold off;

    subplot(2,3,2);
    hold on;
    for ii=1:numDemos
        plot(D(ii).Ys,'k');
    end
    plot(B(2).ymean, 'b');
    plot(xyz_distance(:,2)+B(2).ymean + threshold, 'r');
    plot(-xyz_distance(:,2)+B(2).ymean - threshold, 'r');
    title(sprintf('Set %i - Y', set_num));axis square;grid;
    hold off;

    subplot(2,3,3);
    hold on;
    for ii=1:numDemos
        plot(D(ii).Zs,'k');
    end
    plot(B(2).zmean, 'b');
    plot(xyz_distance(:,3)+B(2).zmean + threshold, 'r');
    plot(-xyz_distance(:,3)+B(2).zmean - threshold, 'r');
    title(sprintf('Set %i - Z', set_num));axis square;grid;
    hold off;

    subplot(2,3,4);
    hold on;
    for ii=1:numDemos
        plot(D(ii).Xs-B(2).xmean,'k');
    end
    plot(xyz_distance(:,1) + threshold, 'b');
    title('Bound X');axis square;grid;
    hold off;
    
    subplot(2,3,5);
    hold on;
    for ii=1:numDemos
        plot(D(ii).Ys-B(2).xmean,'k');
    end
    plot(xyz_distance(:,2) + threshold, 'b');
    title('Bound Y');axis square;grid;
    hold off;
    
    subplot(2,3,6);
    hold on;
    for ii=1:numDemos
        plot(D(ii).Zs-B(2).xmean,'k');
    end
    plot(xyz_distance(:,1) + threshold, 'b');
    title('Bound Z');axis square;grid;
    hold off;
    
    clear d d1 d2 points distance ii threshold R R1 R2
    
    %% calculate and plot canal surfaces
    t = linspace(0,1,nPoints).';
    idx1 = 100;
    idx2 = 100;
    tt = t(idx1:end-idx2);
    xx2 = B(2).xmean(idx1:end-idx2);
    yy2 = B(2).ymean(idx1:end-idx2);
    zz2 = B(2).zmean(idx1:end-idx2);
    RRouter = Router(idx1:end-idx2);
    
    [~,canalouter,~,~,~,~,~] = ...
        canalSurface(tt,xx2,yy2,zz2,RRouter,plotIt,filterIt);

    figure;hold on
    for kk = 1:size(canalouter,3)
        C = canalouter(:,:,kk);
            plot3(C(1,:),C(2,:),C(3,:),'k');
    end
    plot3(xx2,yy2,zz2,'b','linewidth',2);

    for ii=1:numDemos
        plot3(D(ii).Xs,D(ii).Ys,D(ii).Zs,'r','linewidth',2);
    end
    title(sprintf('Set %i - Canal Surface', set_num));
    axis equal
    view(canalview);
    hold off;

end