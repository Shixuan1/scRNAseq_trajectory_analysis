% This script identifies the pseudotime developmental trajectories of hematopoiesis (myeloid lineages) using 10x scRNAseq data of
% bone and bone marrow cells. Multiple branches were detected in this example, including a neutrophil branch, a monocyte/macrophage 
% branch, and a erythroid branch. 
%% load input scRNAseq data
% Change the directory to the current folder to access custom functions (e.g., find_trajectory.m, sort_trajectory.m)
load("data_example2_hematopoiesis.mat"); % 10x scRNAseq data of mouse lemur testis filetered to only include high quality sperm and sperm-progenitor cells. 
% data is organized in a structure, including the following fields
% data.celltypes - a table listing all the cell types involved in the analysis
% data.cells - a table listing all the single cells involved in the analysis
% data.genes - a table listing all the sequenced gene transcripts (including genes not transcribed in the analyzed cells)
% data.mat_raw - a cell by gene matrix listing the transcript count (# of UMI) of each gene in each cell - sparse matrix format is used to save memory

%% data pre-processing
data.mat_libnorm = data.mat_raw./sum(data.mat_raw,2)*1e4; % library size normalization (to 10^4 total transcirpts per cell)
data.mat_lognorm = sparse(log(data.mat_libnorm+1)); % log transformation
data.mat_scaled = normalize(data.mat_lognorm,1); % data scaling by gene, so each gene has mean expression of 0 and std  of 1

%% Identify highly variable genes
% calculate mean and dispersion (VMR) for each gene, scale dispersion conditioned to mean and select genes with high scaled dispersion
xx = mean(data.mat_libnorm,1); % mean expression (average before log transformation)
yy = log(var(data.mat_libnorm,0,1)./xx+1)'; % mean-var-ratio
xx = log(xx+1)'; % log transformation of gene expression level
syy = smooth(xx,yy,0.05,'moving'); % moving average 
yystd = (yy-syy).^2;
syystd = sqrt(smooth(xx,yystd,0.05,'moving')); % moving std
yyscaled = (yy-syy)./syystd; % VMR
data.genes.isvariable = yyscaled>0.5; %*** variable gene threshold

figure(); plot(xx, yy,'.'); hold on
plot(xx(yyscaled>0.5),yy(yyscaled>0.5),'.'); 
plot(xx,syy,'.'); 
plot(xx,syy+syystd*0.5,'.'); hold off;
legend('all genes', 'highly variable genes','moving average','moving threshold')
title("Identify highly variable genes"); xlabel('log mean expression'); ylabel('dispersion (logVMR)')

clear xx yy syy yystd syystd yyscaled

%% PCA for dimensionality reduction
[coeff,score,~,~,explained,~] = pca(full(data.mat_scaled(:,data.genes.isvariable)), 'NumComponents',30); %'Algorithm','eig',
data.PCA = struct;
data.PCA.coeff = coeff;
data.PCA.explained = explained;
data.PCA.score = score;

figure; plot(data.PCA.explained(1:30)/sum(data.PCA.explained(1:30))*100,'o');
xlabel('PC'); ylabel('variance explained, relative among top 30'); title('PCA')

%% Run below to review each top 20 PC to determine if the PC should be used for followup umap generation
% Note that PCs driven by one or a few extreme outliers or immediate early genes should be excluded
isnew = false; % use 'true' here when run the anlaysis the first time or want to select a new sect of PCs; Otherwise, use 'false'; 
if isnew
    mat = data.mat_scaled(:,data.genes.isvariable)';
    genesym = data.genes.symbol(data.genes.isvariable);
    selectPC = true(20,1);
    for nPC = 1:20 % for each PC, visualize genes with top coefficients in top 50 most extreme cells
        [~,idxcell] = sort(data.PCA.score(:,nPC), 'descend');
        [~,idxgene] = sort(data.PCA.coeff(:,nPC), 'descend');
        rkgene = [1:15 length(idxgene)-10:length(idxgene)];
        rkcell = [1:50 length(idxcell)-49:length(idxcell)];
        figure(100); imagesc(mat(idxgene(rkgene),idxcell(rkcell))); colorbar; caxis([-5 5]);
        xlabel('top cells'); ylabel('top genes'); title(["PC"+num2str(nPC); "Press 'y/n' to continue"; "y: use in UMAP, n: removed in UMAP"]);
        yticks(1:length(rkgene)); yticklabels(genesym(idxgene(rkgene)));
        str = input("PC"+num2str(nPC)+", yes/no/break [y/n/x]: ",'s');
        if str=='n'
            selectPC(nPC) = false;
        elseif str=='x'
            disp("Stop considering nPC>"+int2str(nPC));
            selectPC(nPC:end) = false; break;
        elseif str~='y'
            error("Error in manual input! Only 'y/n' is accepted!");
        end
    end
    selectPC = find(selectPC);
else
    selectPC = [1 2 3 4 8 11 13]; % selected PCs
end
data.PCA.selectPC = selectPC(:);    

clear coeff score explained mat idxcell idxgene rkgene rkcell genesym nPC str selectPC isnew

%% UMAP 2D embedding
%opts = struct; opts.perplexity = 25; opts.K = 10;
[~, umap, ~, ~] = run_umap(data.PCA.score(:, data.PCA.selectPC),'Distance','Euclidean');
data.cells.umap = umap.embedding; 

free_annotations = string(unique(data.celltypes.free_annotation)); % unique cell types
colors = [.6 .5 0; 0 0 1; .5 .5 .5; .8 0 1; 0 .6 0; .8 0 0; 0 .8 1; 0 0 0; 0 .8 1; 1 .3 .5; 0 .3 .15; .9 .9 0; 1 .7 .7; .5 .5 1; .3 1 .1]; % colors used to label each cell type

figure; axis equal % visualize UMAP plot
for kk = 1:length(free_annotations)
    idx = data.cells.free_annotation==free_annotations(kk);
    plot(data.cells.umap(idx,1), data.cells.umap(idx,2),'.','Color',colors(kk,:)); hold on
    text(median(data.cells.umap(idx,1)), median(data.cells.umap(idx,2)),strrep(free_annotations(kk), "_", "/"),'Color',colors(kk,:))
end
hold off; box off; axis equal;
%legend(free_annotations,'Location','northeastoutside');
xlim([min(data.cells.umap(:,1))-2 max(data.cells.umap(:,1))+2])
ylim([min(data.cells.umap(:,2))-2 max(data.cells.umap(:,2))+2])
xlabel('UMAP 1'); ylabel('UMAP 2');

clear kk idx mat umap

% save("out_example2_hematopoiesis_1.mat", 'data', '-v7.3'); % recommended save point

%% Set up trajectory identification in UMAP plot (3 trajectories: monocyte/macrophages; netrophil; erythrocyte)
Nc = size(data.cells,1); 
data.cells.trajectory_id = false(Nc,4); 
idx = contains(string(data.cells.free_annotation),["hematopoietic precursor" "granulocyte monocyte progenitor"]); data.cells.trajectory_id(idx,1) = true;
idx = contains(string(data.cells.free_annotation),["granulocyte monocyte progenitor" "neutrophil"]); data.cells.trajectory_id(idx,2) = true;
idx = contains(string(data.cells.free_annotation),["monocyte" "macrophage"]) & ~contains(string(data.cells.free_annotation),...
    "granulocyte monocyte progenitor"); data.cells.trajectory_id(idx,3) = true; %"granulocyte monocyte progenitor" 
idx = contains(string(data.cells.free_annotation),["megakaryocyte progenitor" "hematopoietic precursor" "erythroid"]); data.cells.trajectory_id(idx,4) = true;
clear idx Nc
%% trajectory 1 (hematopoietic-precursor-to-granulocyte-progenitor)
nt = 1; idxc = data.cells.trajectory_id(:,nt);
xx = data.cells.umap(idxc,1); yy = data.cells.umap(idxc,2);
min_x = min(xx); min_y = min(yy); max_x = max(xx); max_y = max(yy);
xrange = [min_x max_x]+(max_x-min_x)/20*[-1 1]; yrange = [min_y max_y]+(max_y-min_y)/20*[-1 1];
[points, xy, distance, t_a] = find_trajectory(xx, yy, 'method', 'auto', 'idx', data.cells.free_annotation(idxc), 'thrDensity', 5,...
    'resolution', min([diff(xrange) diff(yrange)])/500, 'bandwidth', max([diff(xrange) diff(yrange)])/length(xx)*10);  
% data.cells.trajectory(idxc,nt) = t_a;
% data.cells.trajectoryDistance(idxc,nt) = distance;
data.("umap_trajectory"+nt) = points;
%% trajectory 2 (granulocyte-progenitor-to-neutrophils)
nt = 2; idxc = data.cells.trajectory_id(:,nt);
xx = data.cells.umap(idxc,1); yy = data.cells.umap(idxc,2);
min_x = min(xx); min_y = min(yy); max_x = max(xx); max_y = max(yy);
xrange = [min_x max_x]+(max_x-min_x)/20*[-1 1]; yrange = [min_y max_y]+(max_y-min_y)/20*[-1 1];
[points, xy, distance, t_a] = find_trajectory(xx, yy, 'method', 'auto', 'idx', data.cells.free_annotation(idxc), 'thrDensity', 0.9, 'addPoint', data.umap_trajectory1(end,:), ...
    'resolution', min([diff(xrange) diff(yrange)])/25, 'bandwidth', max([diff(xrange) diff(yrange)])/length(xx)*100);  % call find_trajectory function (choose 'auto' or 'manual')
% data.cells.trajectory(idxc,nt) = t_a; 
% data.cells.trajectoryDistance(idxc,nt) = distance;
data.("umap_trajectory"+nt) = points;
%% trajectory 3 (granulocyte-progenitor-to-monocytes/macrophages)
nt = 3; idxc = data.cells.trajectory_id(:,nt);
xx = data.cells.umap(idxc,1); yy = data.cells.umap(idxc,2);
min_x = min(xx); min_y = min(yy); max_x = max(xx); max_y = max(yy);
xrange = [min_x max_x]+(max_x-min_x)/20*[-1 1]; yrange = [min_y max_y]+(max_y-min_y)/20*[-1 1];
[points, xy, distance, t_a] = find_trajectory(xx, yy, 'method', 'auto', 'idx', data.cells.free_annotation(idxc), 'thrDensity', 1.2, 'addPoint', data.umap_trajectory1(end,:), ...
    'resolution', min([diff(xrange) diff(yrange)])/50, 'bandwidth', max([diff(xrange) diff(yrange)])/length(xx)*50);
% data.cells.trajectory(idxc,nt) = t_a; 
% data.cells.trajectoryDistance(idxc,nt) = distance;
data.("umap_trajectory"+nt) = points;
%% trajectory4 (erythrocyte lineage)
nt = 4; idxc = data.cells.trajectory_id(:,nt);
xx = data.cells.umap(idxc,1); yy = data.cells.umap(idxc,2);
min_x = min(xx); min_y = min(yy); max_x = max(xx); max_y = max(yy);
xrange = [min_x max_x]+(max_x-min_x)/20*[-1 1]; yrange = [min_y max_y]+(max_y-min_y)/20*[-1 1];
[points, xy, distance, t_a] = find_trajectory(xx, yy, 'method', 'auto', 'idx', data.cells.free_annotation(idxc), 'thrDensity', 1.2, 'addPoint', data.umap_trajectory1(end,:), ...
    'resolution', min([diff(xrange) diff(yrange)])/50, 'bandwidth', max([diff(xrange) diff(yrange)])/length(xx)*50);  
% data.cells.trajectory(idxc,nt) = t_a; 
% data.cells.trajectory(data.cells.trajectory(:,nt)==0 & contains(string(data.cells.free_annotation),"erythroid")) = ...
%     .5*nanmin(data.cells.trajectory(data.cells.trajectory(:,nt)~=0, nt));
% data.cells.trajectoryDistance(idxc,nt) = distance;
data.("umap_trajectory"+nt) = points;
clear t_a points xy distance xx yy min_x min_y max_x max_y idxc xrange yrange
%% assign trajectory
Nc = size(data.cells,1); 
% match the start of branch 2 and 3 with the end of branch 1
data.umap_trajectory2(1,:) = data.umap_trajectory1(end,:);
data.umap_trajectory3(1,:) = data.umap_trajectory1(end,:);
% calculate trajectory coordinates & distances
data.cells.trajectory = nan(Nc, 6); % trajectory 1-4, 1/2, 1/3
data.cells.trajectoryDistance = nan(Nc,4); 

% recalculate trajectory coordinates, distance for all datapoints 
for nt = 1:4
    [xy,distance,t_a] = distance2curve(data.("umap_trajectory"+nt), data.cells.umap); %[xy,distance,t_a]
    data.cells.trajectoryDistance(:,nt) = distance;
    data.cells.trajectory(:,nt) = t_a;
end
nt = 4; % let MGP be before erythroid
data.cells.trajectory(data.cells.trajectory(:,nt)==0 & contains(string(data.cells.free_annotation),"erythroid")) = ...
    .5*nanmin(data.cells.trajectory(data.cells.trajectory(:,nt)~=0, nt));
nt = 5; % HSC-GPC-neutrophils
[xy,distance,t_a] = distance2curve([data.umap_trajectory1; data.umap_trajectory2], data.cells.umap); 
data.cells.trajectory(:,nt) = t_a;
nt = 6; % HSC-GPC-monocytes
[xy,distance,t_a] = distance2curve([data.umap_trajectory1; data.umap_trajectory3], data.cells.umap); 
data.cells.trajectory(:,nt) = t_a;

[~, data.cells.trajectoryID] = min(data.cells.trajectoryDistance,[],2);
data.cells.istrajOutlier = false(Nc,1);
data.cells.istrajOutlier(~data.cells.trajectory_id((1:Nc)'+Nc*(data.cells.trajectoryID-1))) = true;
data.cells.trajectoryID(data.cells.istrajOutlier) = nan;
clear xy distance t_a nt
% outlierMethod = "moving"; % detect additional outliers
% data.cells.istrajOutlier = false(Nc,1); 
% for nt =1:4
%     switch outlierMethod
%         case "moving"% detect outliers (set moving threshold along trajectory)
%             idxc = find(data.cells.trajectoryID==nt);
%             xx = data.cells.trajectory(idxc,nt);
%             yy = data.cells.trajectoryDistance(idxc,nt);
%             syy = smooth(xx,yy,0.1,'moving'); %smooth(xx,yy,0.1,'rlowess'); %
%             yystd = (yy-syy).^2;
%             syystd = sqrt(smooth(xx,yystd,0.1,'moving')); %sqrt(smooth(xx,yystd,0.1,'rlowess')); %
%             idx = yy>syy+syystd*norminv(0.9999); 
%             figure; plot(xx, yy,'.'); hold on; plot(xx,syy,'.'); plot(xx,syy+syystd*norminv(0.9999),'.'); plot(xx(idx),yy(idx),'o')
%             idx = idxc(idx);
%         case "fixed" % one fixed threshold
%             idx = data.cells.trajectoryDistance>median(data.cells.trajectoryDistance)+mad(data.cells.trajectoryDistance)*1.4826*norminv(0.999);
%     end
%     data.cells.istrajOutlier(idx) = true; pause
% end
% clear nt xx yy yystd syy syystd idxc
%% visualize umap with all trajectory branches
free_annotations = unique(data.celltypes.free_annotation);
figure; 
for kk = 1:length(free_annotations)
    idx = data.cells.free_annotation==free_annotations(kk);%idx = data.cells.free_annotation==unicelltypes(kk);
    plot(data.cells.umap(idx,1), data.cells.umap(idx,2),'.', 'MarkerSize', 10, 'Color', colors(kk,:)); hold on
    text(median(data.cells.umap(idx,1)), median(data.cells.umap(idx,2)), strrep(string(free_annotations(kk)),"_","/"), 'Color', colors(kk,:), 'FontWeight', 'Bold')
end
for nt = 1:4
    [xy,~,~] = distance2curve(data.("umap_trajectory"+nt), data.cells.umap); %[xy,distance,t_a]
    idx = data.cells.trajectoryID==nt & ~data.cells.istrajOutlier;
    plot([data.cells.umap(idx,1)'; xy(idx,1)'],[data.cells.umap(idx,2)'; xy(idx,2)'], 'Color',[1 1 1]*.7, 'HandleVisibility','off');
    plot(data.("umap_trajectory"+nt)(:,1), data.("umap_trajectory"+nt)(:,2), 'k', 'LineWidth',1.5); %fnplt( trajectory,'k'); %
end
    
xlabel('UMAP 1'); ylabel('UMAP 2'); box off; axis equal;
xlim([min(data.cells.umap(:,1))-2 max(data.cells.umap(:,1))+2])
ylim([min(data.cells.umap(:,2))-2 max(data.cells.umap(:,2))+2])

clear kk nt xy idx

%save("out_example2_hematopoiesis_2.mat", 'data', '-v7.3'); % recommended save point

%% calculate correlation of gene expression with trajectory pattern 
% note that this step takes very long time
% Analyze patterns over 3 trajectories: 1) branch 1+ branch 2; 2) branch 1+ branch 3; 3) branch 4; 
npoint = 20; % determines the resolution (1/npoint)
ngenes = size(data.genes,1); %sum(isgene); %size(statsCellType.genes,1);
data.genes.umap_corr = nan(ngenes,3); 
data.genes.umap_fdr = nan(ngenes,3); 
data.genes.umap_peak = nan(ngenes,3); 
nts = [5 6 4]; trajectoryIDs = {[1 2], [1 3], [4]};
for nn = 1:3
    mat_rho = nan(npoint, ngenes);
    mat_pval = nan(npoint, ngenes);
    tpoint = 0:1/npoint:1-1/npoint;
    idxc = ismember(data.cells.trajectoryID, trajectoryIDs{nn}) & ~data.cells.istrajOutlier;
    for kk = 1:length(tpoint)
        [rho, pval] = corr(abs(data.cells.trajectory(idxc,nts(nn))-tpoint(kk)), data.mat_lognorm(idxc,:),'type','Spearman','rows','pairwise');
        mat_rho(kk,:) = rho; mat_pval(kk,:) = pval;
    end
    fdr = mat_pval;% fdr = mafdr(mat_pval(:)); fdr = reshape(fdr, size(mat_pval));
    [~, idx1] = min(fdr,[],1,'linear'); [~, idx2] = min(fdr,[],1');
    rho = mat_rho(idx1)'; fdr =  fdr(idx1)';
    data.genes.umap_corr(:,nn) = rho;
    data.genes.umap_fdr(:,nn) = fdr*(sum(data.genes.isvariable)*npoint);
    data.genes.umap_peak(:,nn) = tpoint(idx2)'; % corresponding t point
end
clear idx1 idx2 mat_rho mat_pval pval kk tpoint fdr rho idxg ia

%save("out_example2_hematopoiesis_3.mat", 'data', '-v7.3'); % recommended save point
%% visualize top genes for each pattern
free_annotations = unique(data.celltypes.free_annotation);

Npattern = 10; thrPval = [1e-300 1e-25 1e-50]; % threshold here is arbitrary and only set to narrow down the number of significant genes
ngenes = size(data.genes,1); 
trajectoryIDs = {[1 2] [1 3] [4]}; nts = [5 6 4];
Nc = size(data.cells,1); data.cells.trajectory_idx = nan(Nc,3); data.genes.trajectory_pattern = nan(ngenes,3);
tabGenes = cell(3,1); tabPattern = cell(3,1);
for mm = 1:3
    % 1. tabGenes
    indc = find(~data.cells.istrajOutlier & ismember(data.cells.trajectoryID,trajectoryIDs{mm}));
    [~, idxc] = sort(data.cells.trajectory(indc, nts(mm)),'ascend');
    
    data.cells.trajectory_idx(indc(idxc),mm) = (1:length(indc))';
    idxg = find(data.genes.umap_fdr(:,mm)<thrPval(mm) & data.genes.isvariable); Ngenes = length(idxg);
    mat = full(data.mat_lognorm(indc(idxc),idxg)');
    smat = nan(size(mat));
    for nn = 1:Ngenes
        smat(nn,:) = smooth(data.cells.trajectory(indc(idxc),nts(mm)), mat(nn,:), 0.1, 'moving'); % 'rlowess'
    end
    idxp = kmeans(smat,Npattern,'Distance','correlation');
    tabGenes{mm} = table; 
    tabGenes{mm}.idxg = idxg; tabGenes{mm}.idxp = idxp; 
    tabGenes{mm}.symbol = data.genes.symbol(tabGenes{mm}.idxg); tabGenes{mm}.ortholog_human = data.genes.ortholog_human(tabGenes{mm}.idxg);
    tabGenes{mm}.mean = -mean(data.mat_lognorm(indc,idxg),1)';
    tabGenes{mm}.fdr = data.genes.umap_fdr(tabGenes{mm}.idxg, mm); tabGenes{mm} = sortrows(tabGenes{mm}, {'idxp', 'fdr', 'mean'});
    tabGenes{mm}.rank = nan(Ngenes,1);
    for nn = 1:Npattern
        idx = tabGenes{mm}.idxp==nn;
        tabGenes{mm}.rank(idx) = (1:sum(idx))';
    end
    figure; imagesc((data.mat_lognorm(indc(idxc),tabGenes{mm}.idxg)./max(data.mat_lognorm(indc(idxc),tabGenes{mm}.idxg),[],2))')
    
    data.genes.trajectory_pattern(tabGenes{mm}.idxg,mm) = tabGenes{mm}.idxp; 
    clear mat smat idxg idx idxp nn
    % 2. tabPattern
    tabPattern{mm} = table;
    tabPattern{mm}.ID = (1:Npattern)';
    tabPattern{mm}.Ngenes = nan(Npattern,1);
    tabPattern{mm}.isReal = false(Npattern,1);
    tabPattern{mm}.peak = nan(Npattern,1);
    tabPattern{mm}.genes = cell(Npattern,1);
    tabGenes{mm}.peak = nan(Ngenes,1); tabGenes{mm}.isReal = false(Ngenes,1);
    for nn = 1:Npattern
        idxg = find(data.genes.trajectory_pattern(:,mm)==nn);
        tabPattern{mm}.Ngenes(nn) = length(idxg);
        tab = table; tab.idxg = idxg; tab.symbol = data.genes.symbol(idxg); tab.fdr = data.genes.umap_fdr(idxg,mm); 
        tab.mean = -mean(data.mat_lognorm(indc(idxc),tab.idxg),1)'; tab = sortrows(tab, {'fdr','mean'});
        pks = median(sum((data.mat_lognorm(indc(idxc),tab.idxg)').*(data.cells.trajectory(indc(idxc),nts(mm))'),2)./sum((data.mat_lognorm(indc(idxc),tab.idxg)'),2)); %mean(data.cells.trajectory.*(data.mat_lognorm(:,idxg)>0.8*max(data.mat_lognorm(:,idxg),[],2)),2);
        tabPattern{mm}.peak(nn) = median(pks); tabGenes{mm}.peak(tabGenes{mm}.idxp==nn) = tabPattern{mm}.peak(nn);
        tabPattern{mm}.genes{nn} = tab;
        figure(40); imagesc(data.mat_lognorm(indc(idxc),tab.idxg)');
        yticks(1:size(tab,1)); yticklabels(tab.symbol);
        str = input("Pattern "+num2str(nn)+", yes/no/break [y/n/x]: ",'s');
        if str=='y'
            tabPattern{mm}.isReal(nn) = true; tabGenes{mm}.isReal(tabGenes{mm}.idxp==nn) = true;
        elseif str~='n'
            error("Error in manual input! Only 'y/n' is accepted!");
        end
    end
    tabPattern{mm} = sortrows(tabPattern{mm}, 'peak', 'ascend');
    tabGenes{mm} = sortrows(tabGenes{mm},{'peak','fdr','mean'});
    
    clear nn pks str tab
end
%% Plot cell type marker genes & top differentially expressed genes along trajectory
Markers = {["CD34" "KIT" "FLT3" "MYB"  ... % hematopoietic precursor
    "STMN1" "MCM5" "TOP2A" "MKI67" ... % proliferation
    "MS4A3" "MPO" "ELANE" "AZU1" "CAMP"  "LTF" "MMP9" "CSF3R" "ALPL"... % neutrophil 
    "S100A8"   ], ... % neutrophil & monocytes 
    ["CD34" "KIT" "FLT3" "MYB"  ... % hematopoietic precursor
    "STMN1" "MCM5" "TOP2A" "MKI67" ... % proliferation
    "S100A8" ... % neutrophil & monocytes
    "IRF8" ... % macrophage & monocytes
    "S100A10" "CSF1R" "VCAN" ... %monocytes
    "MSR1" "C1QA" "APOE" "MRC1" "MARCO"], ... %macrophage
    ["KIT" "MYB" "CD34" "CD38"... % hematopoietic precursor cell
    "STMN1" "MCM5" "TOP2A" "MKI67" ... % proliferation
    "AHSP" "ALAS2"...%erythrocyte lineage
    "GYPA" "XPO7"]}; % mature erythroid

colorMap = [linspace(.8, 1, 100)' linspace(.8, 0, 100)' linspace(.8, 0, 100)';];
trajectoryName = ["precursor to neutrophil" "precursor to monocytes/macrophages" "erythroid"];
for mm = 1:3
    %3. plot
    % plot top 3 detected genes from each cluster
    Nselect = 3; idx = tabGenes{mm}.rank<=Nselect & tabGenes{mm}.isReal;
    tab = table; idxg = tabGenes{mm}.idxg(idx); tab.idxg = idxg;
    tab.pattern = tabGenes{mm}.idxp(idx);
    tab.symbol = tabGenes{mm}.symbol(idx);
    
    indc = find(~data.cells.istrajOutlier & ismember(data.cells.trajectoryID,trajectoryIDs{mm}));
    [~, idxc] = sort(data.cells.trajectory(indc, nts(mm)),'ascend');
    inds = find(contains(string(data.cells.free_annotation(indc(idxc))), "neutrophil")); 
    indS = true(length(idxc),1); indS(inds) = false; indS(inds(1:10:end)) = true; % subsample for neutrophils
    
    mat = data.mat_lognorm(indc(idxc(indS)),idxg)';
    figure; imagesc(mat./max(mat,[],2));
    hold on; yy = find(diff(tab.pattern)~=0); plot([0; Nc], [1; 1]*(yy(:)'+.5), 'k-')
    yticks(1:size(tab,1)); yticklabels(strrep(tab.symbol, "_","\_")); colormap(colorMap); caxis([0 1])
    h = colorbar('Location','eastoutside','Ticks',0:.5:1); ylabel(h,["relative ln(CP10K+1)"]) %,'FontWeight','bold'
    h.Title.String = "gene expression"; title(["Heatmap of top differentially expressed genes" "trajectory "+mm+": "+trajectoryName(mm)])
    xlabel("cells (ordred by trajectory coordinate)"); ylabel("genes (ordred by cluster)");
    
    %plot markers
    markers = Markers{mm};
    [Lia, Locb] = ismember(markers, data.genes.symbol);
    if ~all(Lia); error(join(markers(~Lia))); end
    mat = data.mat_lognorm(indc(idxc(indS)),Locb)';
    figure; imagesc(mat./max(mat,[],2));
    hold on; %yy = find(diff(tab.pattern)~=0); plot([0; Nc], [1; 1]*(yy(:)'+.5), 'k-')
    yticks(1:length(markers)); yticklabels(strrep(markers, "_","\_")); colormap(colorMap); caxis([0 1])
    h = colorbar('Location','eastoutside','Ticks',0:.5:1); ylabel(h,["relative ln(CP10K+1)"]) %,'FontWeight','bold'
    h.Title.String = "gene expression"; title(["Heatmap of marker genes" "trajectory "+mm+": "+trajectoryName(mm)])
    xlabel("cells (ordred by trajectory coordinate)"); ylabel("genes (ordred by cluster)");
    
    % plot cell type bar
    celltypebar = nan(100, length(indc(indS)), 3);
    for nn = 1:length(free_annotations)
        idx = data.cells.free_annotation(indc(idxc(indS)))==free_annotations(nn);
        celltypebar(:,idx,1) = colors(nn,1);
        celltypebar(:,idx,2) = colors(nn,2);
        celltypebar(:,idx,3) = colors(nn,3);
    end
    figure; imshow(celltypebar); title("cell type bar: trajectory "+mm+": "+trajectoryName(mm))
    
    clear xx yy temp nn pks str h idx mat
end


% save output
save("out_example2_hematopoiesis.mat", 'data', 'tabGenes','tabPattern', '-v7.3'); % recommended save point
