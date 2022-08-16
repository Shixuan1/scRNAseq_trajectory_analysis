function points = sort_trajectory(fd_trajectory)
figure; 
thr_branch_length = max(size(fd_trajectory))/10;
fd_trajectory(imfilter(double(fd_trajectory),[1 1 1; 1 0 1; 1 1 1]).*fd_trajectory>=3) = 0;
fd_trajectory = bwlabel(fd_trajectory);
n_curve = max(fd_trajectory(:));
trajectory = table; trajectory.idx = cell(n_curve, 1);
for mm=1:n_curve
    [y_idx, x_idx] = find(fd_trajectory==mm);
    npixel = length(x_idx);
    trajectory.idx{mm} = table;
    trajectory.idx{mm}.y_idx = y_idx;
    trajectory.idx{mm}.x_idx = x_idx;
    trajectory.idx{mm}.id = (1:npixel)';
    trajectory.idx{mm}.type = zeros(npixel,1); %1: end; 2: regular; >2: branching point
    trajectory.idx{mm}.nbranch = zeros(npixel,1); %1: end; 2: regular; >2: branching point
    trajectory.idx{mm}.id_pre = zeros(npixel,1);
    trajectory.idx{mm}.id_post = zeros(npixel,1);
    trajectory.idx{mm}.rank = zeros(npixel,1);
    if length(x_idx)==1
        trajectory.idx{mm}.rank = 1; trajectory.idx{mm}.type = 1;
        continue;
    end
    mat_dist = squareform(pdist([y_idx x_idx]))<=sqrt(2);
    mat_dist(1:npixel+1:end)  = false;
    trajectory.idx{mm}.type = sum(mat_dist,2);
    id_ends = [find(trajectory.idx{mm}.type==1); find(trajectory.idx{mm}.type>2)];
    for kk = 1:length(id_ends)
        if trajectory.idx{mm}.nbranch(id_ends(kk)) == 0
            trajectory.idx{mm}.nbranch(id_ends(kk)) = kk;
        else
            continue;
        end
        id = id_ends(kk); id_pre = 0; count = 1;
        while length(id)==1
            if id_pre~=0
                trajectory.idx{mm}.id_post(id_pre) = id;
            end
            trajectory.idx{mm}.rank(id) = count;
            id_post = find(mat_dist(id,:)'&trajectory.idx{mm}.nbranch==0);
            id_post = id_post(id_post~=id_pre);
            trajectory.idx{mm}.id_pre(id_post) = id;
            trajectory.idx{mm}.nbranch(id) = kk;
            id_pre = id; id = id_post; count = count+1;
            if trajectory.idx{mm}.type(id)>2%||length(id)==1
                trajectory.idx{mm}.rank(id) = count;
                trajectory.idx{mm}.nbranch(id) = kk;
                break
            end
        end
        
    end
    tabcount = groupcounts(trajectory.idx{mm},'nbranch');
    %---------------------VISUALIZAION---------------------------------
   trajectory.idx{mm} = sortrows(trajectory.idx{mm},'rank','ascend');
   %trajectory.idx{mm} = sortrows(trajectory.idx{mm},'nbranch','ascend');
   hold on
   for kk = 1:size(tabcount,1)
       nbranch = tabcount.nbranch(kk);
       if tabcount.GroupCount(kk)<thr_branch_length
           plot(trajectory.idx{mm}.x_idx(trajectory.idx{mm}.nbranch==nbranch),trajectory.idx{mm}.y_idx(trajectory.idx{mm}.nbranch==nbranch),'o-');
       else
           plot(trajectory.idx{mm}.x_idx(trajectory.idx{mm}.nbranch==nbranch),trajectory.idx{mm}.y_idx(trajectory.idx{mm}.nbranch==nbranch),'o-','MarkerFaceColor',[.5 .5 .5]);
       end
   end
   legend; title('filter trajectory'); axis equal; xlim([0 size(fd_trajectory,2)]); ylim([0 size(fd_trajectory,1)])
   trajectory.idx{mm} = sortrows(trajectory.idx{mm},'id','ascend');
   %--------------------------------------------------------------------
%     % remove short branches
%      for kk = 1:size(tabcount,1)
%         if tabcount.GroupCount(kk)<thr_branch_length
%             trajectory.idx{mm}(trajectory.idx{mm}.nbranch== tabcount.nbranch(kk),:) = [];
%         end
%      end   
end

tab = table;
for mm=1:n_curve
    n_trajectory = ones(size(trajectory.idx{mm},1),1)*mm;
    tab = [tab; [trajectory.idx{mm} table(n_trajectory)]];
end
tab = sortrows(tab,'rank','ascend');
uni = unique(tab(:,{'nbranch','n_trajectory'}),'rows');
figure; 
for mm = 1:size(uni,1)
    idx = ismember(tab(:,{'nbranch','n_trajectory'}),uni(mm,:));
    plot(tab.x_idx(idx),tab.y_idx(idx), '.-'); hold on;
end
title('\color{red}Click to select starting point.'); axis equal; xlim([0 size(fd_trajectory,2)]); ylim([0 size(fd_trajectory,1)])

tab.outputrank = nan(size(tab,1),1); tab.segment = nan(size(tab,1),1);
[inputx, inputy] = ginput(1);
dist = sum(abs(([tab.x_idx, tab.y_idx]-[inputx inputy])).^2,2);
idx_pixel = dist == min(dist);
%idx_curve = tab.nbranch==tab.nbranch(idx_pixel)&tab.n_trajectory==tab.n_trajectory(idx_pixel);
count_segment = 1;
tab.segment(tab.nbranch==tab.nbranch(idx_pixel)&tab.n_trajectory==tab.n_trajectory(idx_pixel))=count_segment;
if diff(dist(tab.type~=2&tab.segment==count_segment))<0
    tab.outputrank(tab.segment==count_segment) = sum(tab.segment==count_segment):-1:1;
else
    tab.outputrank(tab.segment==count_segment) = 1:sum(tab.segment==count_segment);
end
tab = sortrows(tab, 'outputrank','ascend');

hold on; plot(tab.x_idx(tab.segment==count_segment),tab.y_idx(tab.segment==count_segment), 'o--'); 
scatter(tab.x_idx(tab.segment==1&tab.outputrank==1),tab.y_idx(tab.segment==1&tab.outputrank==1), 30, 'c','filled');
str = input("Choose next segment of the trajectory? [y/n]: ",'s');
while strcmp(str,'y')
    count_segment = count_segment+1;
    current_rank = max(tab.outputrank);
    title('\color{red}Click to select next segment.'); axis equal; xlim([0 size(fd_trajectory,2)]); ylim([0 size(fd_trajectory,1)])
    [inputx, inputy] = ginput(1);
    dist = sum(abs(([tab.x_idx, tab.y_idx]-[inputx inputy])).^2,2);
    idx_pixel = dist == min(dist);
    %idx_curve = tab.nbranch==tab.nbranch(idx_pixel)&tab.n_trajectory==tab.n_trajectory(idx_pixel);
    tab.segment(tab.nbranch==tab.nbranch(idx_pixel)&tab.n_trajectory==tab.n_trajectory(idx_pixel))=count_segment;
    if ~all(isnan(tab.outputrank(tab.segment==count_segment)))
        warning('Semgent already chosen!'); return;
    end
    dist = sum(abs(([tab.x_idx(tab.segment==count_segment&tab.type~=2), tab.y_idx(tab.segment==count_segment&tab.type~=2)]-...
        [tab.x_idx(tab.outputrank==current_rank) tab.y_idx(tab.outputrank==current_rank)])).^2,2);
    if diff(dist)<0
        tab.outputrank(tab.segment==count_segment) = current_rank+(sum(tab.segment==count_segment):-1:1);
    else
        tab.outputrank(tab.segment==count_segment) = current_rank+(1:sum(tab.segment==count_segment));
    end
    tab = sortrows(tab, 'outputrank','ascend');
    plot(tab.x_idx(tab.outputrank>=current_rank),tab.y_idx(tab.outputrank>=current_rank),'o--')
    if all(~isnan(tab.outputrank))
        disp('All segments assigned.'); break;
    end
    str = input("Choose next segment of the trajectory? [y/n]: ",'s');
end
tab = tab(~isnan(tab.outputrank),:); tab = sortrows(tab, 'outputrank','ascend');
points = [];
for mm = 1:max(tab.segment)
    temp = [tab.x_idx(tab.segment==mm) tab.y_idx(tab.segment==mm)];
    points = [points; temp([1:3:end-1 end],:)];
end

