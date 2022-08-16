function [t_coord, xy, distance, t_a] = find_trajectory(xx,yy, varargin) %method, resolution
%  this function receives 2D scatter datapoints and find a trajectory curve
%  approximately corresponding to the density ridge (method == "auto") or
%  it use manual input to construct a trajectory (method == "manual")
% ----------------INPUT------------------------
% xx: input x axis
% yy: input y axis
% method: find trajectory method ('manual' or 'auto')
% ----------------OUTPUT---------------------
% t_coord, 2D coordinates of the trajectory
% xy: 2D coordinates of the input datapoints projected on the trajectory curve (in the 2D original space)
% distance: distance of the input dataspoints to the trajectory curve
% t_a: 1D coordinates of the input datapoint projected on  the trajectory
%[inputx, inputy] = ginput; %ask for mannual input to guide the curve (in the trajectory space)
%% x y range
min_x = min(xx); min_y = min(yy); max_x = max(xx); max_y = max(yy);
xrange = [min_x max_x]+(max_x-min_x)/20*[-1 1]; yrange = [min_y max_y]+(max_y-min_y)/20*[-1 1];

%% parse input
defaultMethod = 'auto'; expectedMethods = {'auto', 'manual'};
defaultResolution = min([diff(xrange) diff(yrange)])/25; %min([diff(xrange) diff(yrange)])/200; 49
defaultBandwidth = max([diff(xrange) diff(yrange)])/length(xx)*100; %max([diff(xrange) diff(yrange)])/length(xx)*20; 50
defaultThrDensity = 0.9; 
defaultAddPoint = [nan nan];
defaultIdx = ones(length(xx),1); %max([diff(xrange) diff(yrange)])/length(xx)*20; 50
p = inputParser;
addParameter(p,'method',defaultMethod, @(x) any(validatestring(x,expectedMethods)));
addParameter(p,'resolution',defaultResolution, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p,'bandwidth',defaultBandwidth, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p,'thrDensity',defaultThrDensity, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x<10));
addParameter(p,'addPoint',defaultAddPoint, @(x) isnumeric(x) && size(x,2) == 2);
addParameter(p,'idx',defaultIdx, @(x) isnumeric(x) || iscategorical(x));
parse(p,varargin{:}); 

resolution = p.Results.resolution;
bandwidth = p.Results.bandwidth;
thrDensity = p.Results.thrDensity;
addPoint = p.Results.addPoint;
uniidx = unique(p.Results.idx);
%% estimate density
ngrid = [round(diff(xrange)/resolution) round(diff(yrange)/resolution)]; % number of grid for x and y, respectively
grid_x = linspace(xrange(1),xrange(2),ngrid(1));
grid_y = linspace(yrange(1),yrange(2),ngrid(2));
% estimate 2D density distribution (tsne space)
[grid_x1,grid_y1] = meshgrid(grid_x, grid_y);
pts = [grid_x1(:) grid_y1(:)];
[fd,~] = ksdensity([xx yy],pts, 'Bandwidth',bandwidth);
fd = reshape(fd, [ngrid(2) ngrid(1)]);

%% obtain trajectory datapoints from  manual input
if strcmp(p.Results.method,'manual')
    figure; 
    contourf(grid_x, grid_y, fd, 'LineColor',[1 1 1]*.8); hold on % imagesc(grid_x, grid_y, fd_trajectory); hold on
    for nn = 1:length(uniidx)
        idxx = p.Results.idx==uniidx(nn);
        plot(xx(idxx),yy(idxx), '.'); %, 'color', [1 1 1]*.6
    end
    contour(grid_x, grid_y, fd, 'LineColor',[1 1 1]*.8); hold on % imagesc(grid_x, grid_y, fd_trajectory); hold on
    if ~any(isnan(addPoint(:)))
        plot(addPoint(:,1), addPoint(:,2), 'ko--', 'MarkerSize', 10, 'LineWidth', 2); %, 'color', [1 1 1]*.6
    end
    hold off; axis equal; xlim(xrange); ylim(yrange);
    cmap = [linspace(1,1,20)' linspace(1,.5,20)' linspace(1,.5,20)']; colormap(cmap);
    title("\color{red}Left click to input 'additional' trajectory points. Press 'enter' to end.");
    [ridge_x, ridge_y] = ginput; %ask for mannual input to guide the curve
    ridge_x = [addPoint(:,1) ridge_x]; ridge_y = [addPoint(:,2) ridge_y];
end
%% obtain datapoints from automatic desntiy estimation
if strcmp(p.Results.method,'auto')
    %fd_x = reshape(xi(:,1), [ngrid(2) ngrid(1)]); fd_y = reshape(xi(:,2), [ngrid(2) ngrid(1)]);
    [fi,xi] = ksdensity(log(fd(fd(:)~=0)));
    [pks,locs,~,~] = findpeaks(fi,xi);
    fd_trajectory = bwskel(fd>thrDensity*exp(locs(pks == max(pks))));
    if ~any(isnan(addPoint(:)))
        fd_trajectory(abs(addPoint(:,2)-grid_y)==min(abs(addPoint(:,2)-grid_y)), abs(addPoint(:,1)-grid_x)==min(abs(addPoint(:,1)-grid_x))) = true;
    end
    [rows, cols] = find(fd_trajectory);
    points = sort_trajectory(fd_trajectory); % arrange skeleton pixel as a continous curve
    ridge_x = grid_x(points(:,1))'; ridge_y = grid_y(points(:,2))';
    % visualize density estimation & trajectory detection
    figure;  subplot(1,2,1)
    contourf(grid_x, grid_y, fd, 'LineColor',[1 1 1]*.8); hold on % imagesc(grid_x, grid_y, fd_trajectory); hold on
    for nn = 1:length(uniidx)
        idxx = p.Results.idx==uniidx(nn);
        plot(xx(idxx),yy(idxx), '.'); %, 'color', [1 1 1]*.6
    end
    plot(grid_x(cols), grid_y(rows), 'k.'); hold off
    axis xy; axis equal; xlim(xrange); ylim(yrange);title('2D density estimation')
    cmap = [linspace(1,1,20)' linspace(1,.3,20)' linspace(1,.3,20)']; colormap(cmap);
    subplot(1,2,2); 
    %plot(xx,yy, '.', 'color', [.5 .5 .5]); hold on;
    for nn = 1:length(uniidx)
        idxx = p.Results.idx==uniidx(nn);
        plot(xx(idxx),yy(idxx), '.'); hold on %, 'color', [1 1 1]*.6
    end
    contourf(grid_x, grid_y, fd, 'LineColor',[1 1 1]*.8); hold on
    plot(ridge_x, ridge_y, 'r.'); axis equal; xlim(xrange); ylim(yrange); title('trajectory')    
end
%% obtain trajectory points, distance, trajectory coordinates, etc.
%f_trajectory = cscvn( [ridge_x([1:2:end-1 end]), ridge_y([1:2:end-1 end])].' );
f_trajectory = cscvn( [ridge_x([1:end-1 end]), ridge_y([1:end-1 end])].' );
t_coord = fnplt(f_trajectory); t_coord = t_coord';
[xy,distance,t_a] = distance2curve(t_coord, [xx yy]); %[xy,distance,t_a]
hold on; plot(t_coord(:,1),t_coord(:,2), 'k-'); %plot([xx'; xy(:,1)'],[yy'; xy(:,2)'],'color',[1 1 1]*.5)

