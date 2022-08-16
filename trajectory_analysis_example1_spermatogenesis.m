% This script identifies the pseudotime developmental trajectory of spermatogenesis using 10x scRNAseq data of sperm and sperm progenitor
% cells, and detect the genes that are differentially expressed along the trajectory
%% load input scRNAseq data
% Change the directory to the current folder to access custom functions (e.g., find_trajectory.m, sort_trajectory.m)
load("data_example1_spermatogenesis"); % 10x scRNAseq data of mouse lemur testis filetered to only include high quality sperm and sperm-progenitor cells. 
% data is organized in a structure, including the following fields
% data.celltypes - a table listing all the cell types involved in the analysis
% data.cells - a table listing all the single cells involved in the analysis
% data.genes - a table listing all the sequenced gene transcripts (including genes not transcribed in the analyzed cells)
% data.mat_raw - a cell by gene matrix listing the transcript count (# of UMI) of each gene in each cell

%% data pre-processing
data.mat_libnorm = data.mat_raw./sum(data.mat_raw,2)*1e4; % library size normalization (to 10^4 total transcirpts per cell)
data.mat_lognorm = log(data.mat_libnorm+1); % log transformation
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
[coeff,score,~,~,explained,~] = pca(data.mat_scaled(:,data.genes.isvariable), 'NumComponents',30); %'Algorithm','eig',
data.PCA = struct;
data.PCA.coeff = coeff;
data.PCA.explained = explained;
data.PCA.score = score;

figure; plot(data.PCA.explained(1:30)/sum(data.PCA.explained(1:30))*100,'o');
xlabel('PC'); ylabel('variance explained, relative among top 30'); title('PCA')

% Run below to review each top 20 PC to determine if the PC should be used for followup umap generation
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
    selectPC =  1:4; % selected PCs
end
data.PCA.selectPC = selectPC(:);    

clear coeff score explained mat idxcell idxgene rkgene rkcell genesym nPC str selectPC isnew

%% UMAP 2D embedding
%opts = struct; opts.perplexity = 25; opts.K = 10;
[~, umap, ~, ~] = run_umap(data.PCA.score(:, data.PCA.selectPC),'Distance','Euclidean');
data.cells.umap = umap.embedding; 

colors = [0 .8 1; 0 0 1; 0 .6 0; .6 .5 0; .8 0 0; .8 0 1; .5 .5 .5]; % colors used to label each cell type

figure; axis equal % visualize UMAP plot
for kk = 1:size(data.celltypes,1)
    idx = data.cells.free_annotation==data.celltypes.free_annotation(kk);
    plot(data.cells.umap(idx,1), data.cells.umap(idx,2),'.','Color',colors(kk,:)); hold on
    text(median(data.cells.umap(idx,1)), median(data.cells.umap(idx,2)),strrep(string(data.celltypes.free_annotation(kk)), "_", "/"),'Color',colors(kk,:))
end
hold off; box off; axis equal;
%legend(free_annotations,'Location','northeastoutside');
xlim([min(data.cells.umap(:,1))-2 max(data.cells.umap(:,1))+2])
ylim([min(data.cells.umap(:,2))-2 max(data.cells.umap(:,2))+2])
xlabel('UMAP 1'); ylabel('UMAP 2');

clear kk idx mat umap

%% Find trajectory in UMAP plot
Nc = size(data.cells,1); data.cells.trajectory = nan(Nc,1); 
data.umap_trajectory = nan(Nc,2); 
data.cells.istrajOutlier = false(Nc,1);

% Call find_trajectory function (choose 'auto' or 'manual')
% Note that in this specific example, we selected 'auto' (automatic detection of trajectory) and default parameters (see parameter choice in function 'find_trajectory')
% Two discountious trajectories have been detected, a single dot at the start of the trajectory and the rest of the trajectory, please follow the
% promted figure to assign the direction of the trajectory and connect the two segments
[points, xy, distance, t_a] = find_trajectory(data.cells.umap(:,1), data.cells.umap(:,2),'method', 'auto');  
data.cells.trajectory(:) = t_a;
data.umap_trajectory = points;
data.cells.trajectory_distance = distance;

% detect outlier datapoints % note that in this specific example, no cell was deemed outlier
outlierMethod = "moving";
switch outlierMethod % two methods were provided to detect outliers
    case "moving"% detect outliers (set moving threshold along trajectory, suitable for data with different dispersion along the trajectory)
        xx = data.cells.trajectory(:);
        yy = data.cells.trajectory_distance;
        syy = smooth(xx,yy,0.25,'moving'); 
        yystd = (yy-syy).^2;
        syystd = sqrt(smooth(xx,yystd,0.25,'moving'));  
        idx = yy>syy+syystd*norminv(0.9999); % change here to relax or tighten the threshold
        figure; plot(xx, yy,'.'); hold on; plot(xx,syy,'.'); plot(xx,syy+syystd*norminv(0.9999),'.'); plot(xx(idx),yy(idx),'o')
        title('Detect outliers'); xlabel('trajectory'); ylabel('distance to the trajectory')
        legend('cells', 'moving average','moving threshold','outliers')
    case "fixed" % one fixed threshold
        thr = median(distance)+mad(distance)*1.4826*norminv(0.999); % change here to relax or tighten the threshold
        idx = distance>thr;
        figure; plot(xx, yy,'.'); hold on; plot([0 1], [1 1]*thr, '.'); plot(xx(idx),yy(idx),'o')
        title('Detect outliers'); xlabel('trajectory'); ylabel('distance to the trajectory')
        legend('cells', 'moving average','moving threshold','outliers')
end
data.cells.istrajOutlier(idx) = true; 

figure; % visualize UMAP plots with trajectory
for kk = 1:size(data.celltypes,1)
    idx = data.cells.free_annotation==data.celltypes.free_annotation(kk);%idx = data.cells.free_annotation==unicelltypes(kk);
    plot(data.cells.umap(idx,1), data.cells.umap(idx,2),'.', 'MarkerSize', 10, 'Color', colors(kk,:)); hold on
    text(median(data.cells.umap(idx,1)), median(data.cells.umap(idx,2)), data.celltypes.free_annotation(kk), 'Color', colors(kk,:), 'FontWeight', 'Bold')
end
if  exist('xy','var')~=1
    [xy,~,~] = distance2curve(data.umap_trajectory, [data.cells.umap(:,1) data.cells.umap(:,2)]); %[xy,distance,t_a]
end
idx = ~data.cells.istrajOutlier;
plot([data.cells.umap(idx,1)'; xy(idx,1)'],[data.cells.umap(idx,2)'; xy(idx,2)'], 'Color',[1 1 1]*.7, 'HandleVisibility','off');
%legendlabel = celltypelabels;
if exist('inputx','var')==1 && exist('inputy','var')==1
    plot(inputx,inputy,'ko'); %legendlabel = [legendlabel; "trajectory mannual input"];
end
plot(data.umap_trajectory(:,1), data.umap_trajectory(:,2), 'k', 'LineWidth',1.5); %fnplt( trajectory,'k'); %
    
xlabel('UMAP 1'); ylabel('UMAP 2'); box off; axis equal;
xlim([min(data.cells.umap(:,1))-2 max(data.cells.umap(:,1))+2])
ylim([min(data.cells.umap(:,2))-2 max(data.cells.umap(:,2))+2])

clear xx yy syy yystd syystd idx thr xy t_a kk idx outlierMethod points distance

% save("output_example1_spermatogenesis_1.mat", 'data', '-v7.3') % recommended point to save the trajectory data

%% Detect genes following a trajectory expression pattern
% this analysis compares each gene's expression pattern along the trajectory with 20 unimodal patterns 
% (with their single peaks uniformaly distributed from the beginning to the end of the trajectory)
npoint = 20; % determines the resolution (1/npoint)
ngenes = size(data.genes,1); 
mat_rho = nan(npoint, ngenes); % to store correlation score
mat_pval = nan(npoint, ngenes); % to store pval
tpoint = 0:1/npoint:1-1/npoint;
for kk = 1:length(tpoint)
    [rho, pval] = corr(abs(data.cells.trajectory(~data.cells.istrajOutlier)-tpoint(kk)), data.mat_lognorm(~data.cells.istrajOutlier,:),'type','Spearman','rows','pairwise');
    mat_rho(kk,:) = rho; mat_pval(kk,:) = pval;
end
fdr = mat_pval;
[~, idx1] = min(fdr,[],1,'linear'); [~, idx2] = min(fdr,[],1');
rho = mat_rho(idx1)'; fdr =  fdr(idx1)';
% the pattern with the best correlation is kept
data.genes.umap_corr = nan(size(data.genes,1),1); data.genes.umap_corr = rho; 
data.genes.umap_fdr = nan(size(data.genes,1),1); data.genes.umap_fdr = fdr*(sum(data.genes.isvariable)*npoint); % multi-testing correction - note that only highly-variable genes were considered
data.genes.umap_peak = nan(size(data.genes,1),1); data.genes.umap_peak = tpoint(idx2)'; % the pattern of the best correlation
data.genes.umap_rk_fdr = nan(size(data.genes,1),1); data.genes.umap_rk_fdr = tiedrank(data.genes.umap_fdr(~isnan(data.genes.umap_peak)),'ascend');

clear idx1 idx2 mat_rho mat_pval pval kk tpoint fdr rho idxg ia npoint ngenes

%% Cluster the significant genes (k-means clustering)
Npattern = 10; % number of k-means cluster
thrPval = 1e-100; % pvalue threshold
Nc = size(data.cells,1); ngenes = size(data.genes,1); 
indc = find(~data.cells.istrajOutlier);
[~, idxc] = sort(data.cells.trajectory(indc),'ascend');
data.cells.trajectory_idx = nan(Nc,1);
data.cells.trajectory_idx(indc(idxc)) = (1:length(indc))';
idxg = find(data.genes.umap_fdr<thrPval & data.genes.isvariable); % select the significant highly variable genes for downstream analysis
Ngenes = length(idxg);
mat = data.mat_lognorm(indc(idxc),idxg)';
smat = nan(size(mat)); % use smoothed expression pattern in clustering to reduce noise
for nn = 1:Ngenes
    smat(nn,:) = smooth(data.cells.trajectory(indc(idxc)), mat(nn,:), 0.1, 'moving'); % 'rlowess'
end
idxp = kmeans(smat,Npattern,'Distance','correlation'); % k-means clustering

tabGenes = table; tabGenes.idxg = idxg; tabGenes.cluster = idxp; 
tabGenes.symbol = data.genes.symbol(tabGenes.idxg); tabGenes.ortholog_human = data.genes.ortholog_human(tabGenes.idxg);
tabGenes.fdr = data.genes.umap_fdr(tabGenes.idxg); tabGenes = sortrows(tabGenes, {'cluster', 'fdr'});
tabGenes.rank = nan(Ngenes,1);
for nn = 1:Npattern
    idx = tabGenes.cluster==nn;
    tabGenes.rank(idx) = (1:sum(idx))';
end
data.genes.trajectory_pattern = nan(ngenes,1);
data.genes.trajectory_pattern(tabGenes.idxg) = tabGenes.cluster; %all(data.genes.Symbol(tabGenes.idxg)==tabGenes.Symbol)

figure; imagesc((data.mat_lognorm(indc(idxc),tabGenes.idxg)./max(data.mat_lognorm(indc(idxc),tabGenes.idxg),[],2))'); colorbar; caxis([0 0.85]); colormap('hot')
title('Heatmap of significant gene expression'); xlabel('cells (ordered by trajectory coordinate)'); ylabel('genes (ordered by cluster)')

clear mat smat idxg idx idxp nn Nc ngenes thrPval

%======== Review each gene cluster===========
tabPattern = table; 
tabPattern.ID = (1:Npattern)'; 
tabPattern.Ngenes = nan(Npattern,1);
tabPattern.isReal = false(Npattern,1);
tabPattern.peak = nan(Npattern,1);
tabPattern.genes = cell(Npattern,1);
tabGenes.peak = nan(Ngenes,1); tabGenes.isReal = false(Ngenes,1);
for nn = 1:Npattern
    idxg = find(data.genes.trajectory_pattern==nn);
    tabPattern.Ngenes(nn) = length(idxg);
    tab = table; tab.idxg = idxg; tab.Symbol = data.genes.symbol(idxg); tab.fdr = data.genes.umap_fdr(idxg); tab = sortrows(tab, 'fdr');
    pks = median(sum((data.mat_lognorm(indc(idxc),tab.idxg)').*(data.cells.trajectory(indc(idxc))'),2)./sum((data.mat_lognorm(indc(idxc),tab.idxg)'),2)); %mean(data.cells.trajectory.*(data.mat_lognorm(:,idxg)>0.8*max(data.mat_lognorm(:,idxg),[],2)),2);
    tabPattern.peak(nn) = median(pks); tabGenes.peak(tabGenes.cluster==nn) = tabPattern.peak(nn);
    tabPattern.genes{nn} = tab;
    figure(40); imagesc(data.mat_lognorm(indc(idxc),tab.idxg)'); 
    title("Heatmap of cluster " +nn+" genes"); xlabel('cells (ordered by trajectory coordinate)'); ylabel('genes (ordered by cluster)')

    yticks(1:size(tab,1)); yticklabels(tab.Symbol);
    % review if the cluster is correctly detected
    str = input("Pattern "+num2str(nn)+", yes/no/break [y/n/x]: ",'s');
    if str=='y'
        tabPattern.isReal(nn) = true; tabGenes.isReal(tabGenes.cluster==nn) = true;
    elseif str~='n'
        error("Error in manual input! Only 'y/n' is accepted!");
    end
end
tabPattern = tabPattern(tabPattern.isReal,:);
tabPattern = sortrows(tabPattern, 'peak', 'ascend'); % order clusters by average peak
tabGenes = sortrows(tabGenes,{'peak','fdr'});

clear nn pks str tab 

% save all analysis output
save("output_example1_spermatogenesis.mat", 'data', 'tabGenes', 'tabPattern', '-v7.3');

% ===========plot top 5 genes from each cluster===============
Nselect = 5; idx = tabGenes.rank<=Nselect & tabGenes.isReal;
tab = table; idxg = tabGenes.idxg(idx); tab.idxg = idxg; 
tab.pattern = tabGenes.cluster(idx);
tab.Symbol = tabGenes.symbol(idx);  

colorMap = [linspace(.8, 1, 100)' linspace(.8, 0, 100)' linspace(.8, 0, 100)';];
mat = data.mat_lognorm(indc(idxc),idxg)';
figure; imagesc(mat./max(mat,[],2));
hold on; yy = find(diff(tab.pattern)~=0); plot([0; Nc], [1; 1]*(yy(:)'+.5), 'k-')
yticks(1:size(tab,1)); yticklabels(strrep(tab.Symbol, "_","\_")); colormap(colorMap); caxis([0 1])
h = colorbar('Location','eastoutside','Ticks',0:.5:1); ylabel(h,["relative ln(CP10K+1)"]) %,'FontWeight','bold'
h.Title.String = "gene expression"; title("Heatmap of top " +Nselect+ " genes from each cluster")
xlabel("cells (ordered by trajectory coordinate)"); ylabel("genes (ordered by cluster)")

%% Plot expression of marker genes along the trajectory
% examples of spermatogenesis markers
markers = ["ID4" "UCHL1" "SALL4" "KIT" "DMRT1" "ZBTB16" "CCNA2" "GFRA1" "TJP3" "SOHLH1" "RHOXF1" "OVOL2" "CRABP1" ...
     "MLH1" ...
     "PIWIL2" "STRA8" "SYCP2" ... 
     "SYCP3" "SYCP1" "SPAG6" "SYCE1" "PTTG1" "HORMAD1" "TBPL1" "CCNA1" "INSL6"  ...
     "DDX4" "ACR" "NME8" "SPO11" "ACRV1" ...
    "TEX101" "TSSK1B" "TXNDC2" "TNP1" "AZIN2" "PRM1" "GAPDHS" "OAZ3" "PGK2"];

[Lia, Locb] = ismember(markers, data.genes.symbol);
if ~all(Lia); error(join(markers(~Lia))); end
mat = data.mat_lognorm(indc(idxc),Locb)';
figure; imagesc(mat./max(mat,[],2));
hold on; %yy = find(diff(tab.pattern)~=0); plot([0; Nc], [1; 1]*(yy(:)'+.5), 'k-')
yticks(1:length(markers)); yticklabels(strrep(markers, "_","\_")); colormap(colorMap); caxis([0 1])
h = colorbar('Location','eastoutside','Ticks',0:.5:1); ylabel(h,["relative ln(CP10K+1)"]) %,'FontWeight','bold'
h.Title.String = "gene expression"; title("Heatmap of marker genes")
xlabel("cells (ordered by trajectory coordinate)"); ylabel("genes (ordered by cluster)")

% plot cell type bar
%colors = [0 .8 1; 0 0 1; 0 .6 0; .6 .5 0; .8 0 0; .8 0 1; .5 .5 .5]; % colors used to label each cell type  %colors = repmat([.6 .5 0; 0 0 1; .5 .5 .5; .8 0 1; 0 .6 0; .8 0 0; 0 .8 1; 0 0 0; 1 .5 0; 1 .2 .2],3,1);

celltypebar = nan(100, length(indc), 3);
for nn = 1:size(data.celltypes,1)
    idx = data.cells.free_annotation(indc(idxc))==data.celltypes.free_annotation(nn);
    celltypebar(:,idx,1) = colors(nn,1);
    celltypebar(:,idx,2) = colors(nn,2);
    celltypebar(:,idx,3) = colors(nn,3);
end
figure; imshow(celltypebar); title("cell type")

clear xx yy temp nn pks str h idx mat