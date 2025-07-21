clear
master_dir = pwd
addpath(master_dir)
addpath([master_dir,'/Util'])
addpath(genpath('/data/wheelock/data1/people/Cindy/PublicRepo/BrBx-HSB_infomap_cleanup'));
rmpath('/data/wheelock/data1/people/Cindy/PublicRepo/BrBx-HSB_infomap_cleanup/tmp');
rmpath('/data/wheelock/data1/people/Cindy/PublicRepo/BrBx-HSB_infomap_cleanup/legacy');
addpath(genpath('/data/wheelock/data1/people/Cindy/PublicRepo/NLA_toolbox_070319/NLA visualizationfunctions'));

load('/data/wheelock/data1/people/Cindy/PublicRepo/GeneralizabilityAreaParcellation/NetworkCommunityAssignments/IM_Tu326_19Networks.mat','IM');
IM_Tu19 = IM;
load('/data/wheelock/data1/people/Cindy/PublicRepo/GeneralizabilityAreaParcellation/NetworkCommunityAssignments/IM_Tu326_12Networks.mat','IM');
IM_Tu12 = IM;

load('IM_Gordon_13nets_333Parcels_Renamed2024.mat','IM');
IM_Gordon = IM;
load('/data/wheelock/data1/parcellations/IM/IM_MyersLabonte_23Networks.mat','IM');
IM_MyersLabonte = IM; IM_MyersLabonte.Nets{end} = 'Unspecified';
load('/data/wheelock/data1/parcellations/IM/IM_Gordon2017vertex','IM');
IM_Gordon2017vertex = IM;

clear IM
% Load parcels
Parcels_Tu_326 = smartload('Parcels_Tu_326.mat');
Parcels_Gordon = smartload('Parcels_Gordon.mat');
Parcels_MyersLabonte283 = smartload('Parcels_Myers_Labonte_202310.mat');

load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;clear MNIl MNIr

namestr1 = 'eLABE_Y2_N113_atleast600frames'
namestr2 = 'eLABE_Y3_N86'
parcelname = 'Tu_326'
IM = IM_Tu19;
Parcels = Parcels_Tu_326;
K = 19
centroid_filename = 'eLABE_Y2_Y3_template_dim326'

[SI,h] = silhouette(template_X,template_y,'correlation');
mSI = mean(SI)
figure('position',[100 100 250 250]);
histogram(SI,-1:0.05:1,'EdgeColor','none','normalization','PDF');
xline(mSI,'Color','r');
xlabel('SI');
ylabel('PDF');
set(gca,'FontSize',12);
print(['./Figure/eLABE_template_initial_SI_',method,'Zdim',num2str(zdim)],'-dpng','-r300');

figure('position',[100 100 800 400]);
xlim([0,K+1]);
ylim([-1,1]);
for inet = 1:K
    v = Violin(SI(template_y==inet),inet,'ViolinColor',IM.cMap(inet,:));
    v.ScatterPlot.SizeData = 10;
end
ylabel('SI');
xlabel('Networks');
xticks(1:max(template_y));
xticklabels(IM.Nets);
xtickangle(35);
set(gca,'FontSize',12);
print(['./Figure/eLABE_template_initial_SI_violinplot_',method,'Zdim',num2str(zdim)],'-dpng','-r300');

centroids_init_corr = NaN(K,nROI);D_init = cell(K,1);sumd_init = NaN(K,1);
for inet = 1:K
    [~,centroids_init_corr(inet,:),sumd_init(inet),D_init{inet}] = kmeans(template_X(template_y==inet,:),1,'distance','correlation','start','plus','Display','off','Replicates',1);
end
%% Plot the intial centroids based on the group-level assignments
cmap = ROY_BIG_BL(100);
colorrange = [-0.2,0.2];
for inet = 1:K
    vals = centroids_init_corr(inet,:)';
    f = figure('position',[100 100 385 275]);
    ax1 = subplot(2,1,1);
    set(ax1,'Position',[0,0.5,0.85,0.5]);
    plot_parcels_by_values(vals,Anat,'lat',Parcels,colorrange,cmap)
    ax2 = subplot(2,1,2);
    set(ax2,'Position',[0,0.05,0.85,0.5]);
    plot_parcels_by_values(vals,Anat,'med',Parcels,colorrange,cmap)
    
    h = axes(f,'visible','off'); % attach colorbar to h
    c = colorbar(h,'Position',[0.88 0.1680 0.022 0.7],'XTick',[0,1],'XTicklabel',colorrange,'FontSize',12);
    colormap(c,cmap);
    exportgraphics(gcf,['./Figure/NetworkTemplate_',num2str(inet),'.png'],'Resolution',300)
    close all
end

%% Plot the distance to initial centroids
N = ceil(sqrt(K));
figure('position',[100 100 800 800]);clear ax
t = tiledlayout(N,N,'TileSpacing','tight');
thre = NaN(1,K);
for inet = 1:K
    ax(inet) = nexttile;
    histogram(1-D_init{inet},'normalization','PDF','EdgeColor',[0.5 0.5 0.5],'FaceColor',IM.cMap(inet,:),'FaceAlpha',1);
    set(gca,'FontSize',12);
    yticks([]);xticks([]);
    xlim([-1,1])
end
linkaxes(ax,'xy');
ylabel(t,'PDF','FontSize',12);
xlabel(t,'Correlation to template','FontSize',12);
print('./Figure/Correlation_to_init_centroids.png','-dpng','-r300');

% Across all networks
figure('position',[100 100 200 200]);
histogram(1-cat(1,D_init{:}),'normalization','PDF','EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',1);
ylabel('PDF');
xlabel('r');
xlim([-1,1])
set(gca,'FontSize',12);
print('./Figure/Distance_to_init_centroids.png','-dpng','-r300');

% Report mean and SD
tmp = atanh(1-cat(1,D_init{:}));
tanh(mean(tmp))
tanh(std(tmp))

%% Save centroids
save([centroid_filename,'_centroids.mat'],'centroids*');

%% K-means with cross-validation to find K


%% Calculate different number of clusters
% N.B. It seems 100 replicates produces more stable results across split
% halves than 10 replicates
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);
if isempty(gcp('nocreate'))
    parpool(10);
end
tic
for k = 1:maxk
    for isplit = 1:Nsplits
        split1_id = false(nROI,Nsess);split1_id(:,split1_sub{isplit})=true;
        split1_id = reshape(split1_id,[],1);
        rng(1);
        [split1_assn,split1_C] = kmeans(template_X(split1_id,:),k,'distance','correlation','start','plus','Display','final','Replicates',100,'options',options);
        [~,split2_hat] = pdist2(split1_C,template_X(~split1_id,:),'correlation','Smallest',1); % Find the nearest centroid
        rng(1);
        [split2_assn] = kmeans(template_X(~split1_id,:),k,'distance','correlation','start','plus','Display','final','Replicates',100,'options',options);
        [ARI(isplit,k)] = adj_rand_index_mod(split2_hat,split2_assn);
        [~,NMI(isplit,k)] = partition_distance(split2_hat',split2_assn);
    end
    toc
end
save([centroid_filename,'_kmeans_stability.mat'],'ARI','NMI')

%% Plot stability
load([centroid_filename,'_kmeans_stability.mat'],'ARI','NMI')

figure;
boxplot(ARI)
ylim([0.5,1]);
yticks([0.5:0.25:1]);
xticks([0:5:maxk]);
xticklabels([0:5:maxk]);
ylabel('ARI');
xlabel('Number of clusters');
grid on
set(gca,'FontSize',12);
% print(['./Figure/',centroid_filename,'_kmeans_ARI'],'-dpng','-r300')

figure;
boxplot(NMI)
ylim([0.5,1]);
yticks([0.5:0.25:1]);
xticks([0:5:maxk]);
xticklabels([0:5:maxk]);
ylabel('NMI');
xlabel('Number of clusters');
grid on
set(gca,'FontSize',12);
% print(['./Figure/',centroid_filename,'_kmeans_NMI'],'-dpng','-r300')

%% Obtain final best k


%% Fit the final k with all data and save the best k centroid
bestk = 23 % based on stability
tic
rng(1);[assn_bestk,C_bestk] = kmeans(template_X,bestk,'distance','correlation','start','plus','Display','final','Replicates',1000,'options',options);
toc

%% if we have the Centroids already
[D_bestk,assn_bestk] = pdist2(C_bestk,template_X,'correlation','Smallest',1);
assn_bestk = assn_bestk';

%% View on the brain
k = bestk
eval('assn = assn_bestk;',num2str(bestk));

assn = reshape(assn,nROI,[]);
prop = calc_network_prop(assn,1:k);

for iNet = 1:k%k
    View_prop_Colors_transparent(prop(:,iNet),Parcels)
    %     print(gcf,fullfile('./Figure/',sprintf('ProportionNetworks%02d',iNet)),'-dpng','-r300');
    %     pause;
    close all;
end

% Get consensus across template subjects
[~,consensus_bestk] = max(prop,[],2);
plot_network_assignment_parcel_key(Parcels,consensus_bestk);
print(['./Figure/consensus_assn_bestk.png'],'-dpng','-r300')
try
    CW = smartload(['./template_matching_results/IM_consensus_eLABE_Y2_Y3_assn_bestk',num2str(bestk),'.mat']);
catch
    switch bestk
        case 23
            %% Assign colors for 23 network solution
            [CW,GenOrder] = assign_networks_by_template(consensus_bestk,IM_Tu19,0.3);
            CW.Nets{1} = 'Language'; CW.cMap(1,:) = IM_Tu19.cMap(10,:);
            CW.Nets{2} = 'Frontoparietal';CW.cMap(2,:) =[153 101 21]/255;% [255 223 0]/255;
            % CW.Nets{4}= 'Primary Visual'; CW.cMap(4,:) =[153 101 21]/255;
            CW.Nets{5}= 'Parietal Memory';
            CW.Nets{6}= 'Ventral Insula'; CW.cMap(6,:) = [0.7 0.8 0.7];
            CW.Nets{7} = 'Premotor/Dorsal Attention II'
            CW.Nets{8} = 'Inferior Orbitofrontal';CW.cMap(8,:) = [240,255,100]/255;%
            CW.Nets{10}= 'Medial Default';
            CW.Nets{11}= 'Context Association';
            CW.Nets{12}= 'Somato-cognitive-action'; CW.cMap(12,:) = [126 0 102]/255;
            CW.Nets{13} = 'Default-Parietal';
            CW.Nets{14} = 'Visual Peripheral';
            CW.Nets{15}= 'Action Mode';
            CW.Nets{17} = 'Visual-Streams';
            CW.Nets{18} = 'Default-Anterolateral';
            CW.Nets{19} = 'Salience II';CW.cMap(19,:) = [132 176 103]/255;
            CW.Nets{20} = 'Temporal Default';CW.cMap(20,:) = [255 0 0]/255;
            CW.Nets{23} = 'Dorsal SM';CW.cMap(23,:) = IM_Tu19.cMap(18,:);
        case 7
            [CW,GenOrder] = assign_networks_by_template(consensus_bestk,IM_Tu12,0.3);
        case 17
            IM23 = load('./template_matching_results/IM_consensus_eLABE_Y2_Y3_assn_bestk23.mat');
            [CW,GenOrder] = assign_networks_by_template(consensus_bestk,IM23,0.3);
            CW.Nets{13} = 'Dorsal SM'; CW.cMap(13,:) = IM_Tu19.cMap(18,:);
    end
end
inet =6
plot_network_assignment_parcel_key(Parcels,double(consensus_bestk==inet),CW.cMap(inet,:));

plot_network_assignment_parcel_key(Parcels,consensus_bestk,CW.cMap);
print(gcf,['./Figure/consensus_assn_bestk.png'],'-dpng','-r300')

% Plot a legend for the Networks
N = length(CW.Nets)
figure('Units','inches','position',[10 10 8 2])%[10 10 5,3]);%[10 10 5 2] %[10 10 6,3])
h = gscatter(ones(1,N),ones(1,N),CW.Nets,CW.cMap,'s',50);
for i = 1:N
    set(h(i),'Color','k','MarkerFaceColor',CW.cMap(i,:));
end
legend(CW.Nets,'interpreter','none','FontSize',10,'location','best','Orientation','horizontal','NumColumns',3);
% legend(CW.Nets,'interpreter','none','FontSize',10,'location','best','Orientation','horizontal','NumColumns',1);
legend('boxoff')
xlim([10,11]);
axis('off')
print('./Figure/networks_Legend','-dpng','-r300');

%% Plot network clusters
inet = 1;
figure;
imagesc(template_X(assn_bestk==inet,:));
colormap(ROY_BIG_BL);caxis([-0.5,0.5]);
axis off
print('./Figure/network_FC_concatenated','-dpng','-r300');

return

%% Write to IM file that can be loaded to plot the networks on the brain
clear IM
IM.name=['IM_consensus_eLABE_Y2_Y3_assn_bestk',num2str(bestk)];
IM.cMap=CW.cMap;
IM.Nets= CW.Nets;
load('/data/wheelock/data1/people/Cindy/BCP/ParcelPlots/Parcels_Tu_326.mat','ROIxyz');
IM.ROIxyz=ROIxyz;
Nroi = length(consensus_bestk);
IM.key=[[1:Nroi]',zeros(Nroi,1)];
[IM.key(:,2),IM.order]=sort(IM_Remove_Naming_Gaps_HSB(consensus_bestk));
IM.ROIxyz=IM.ROIxyz(IM.order,:);
IM=Org_IM_DVLR_HSB(IM);
save(['./template_matching_results/',IM.name,'.mat'],'IM');

% Plot abundance across networks
[~,sort_sub_id] = sort(T.age_yrs(T.in_template==1));

figure;
imagesc(assn(IM.order,sort_sub_id));
hl = xline(40.5); % where eLABE Y2/Y3 splits
hl.LineWidth = 2;
colormap(IM.cMap);
xticks([]); yticks([]);
xlabel('Session (sorted young to old)');
ylabel('Parcel');
set(gca,'FontSize',12);
print('./Figure/assn_across_subjects.png','-dpng','-r300');

%% Plot individualized networks
k = bestk
assn = assn_bestk;
cmap = CW.cMap

assn = reshape(assn,nROI,[]);

icol = 28 %14/42 28/66
plot_network_assignment_parcel_key(Parcels, assn(:,icol),cmap)
print(['./Figure/brain_plot_session',num2str(icol),'.png'],'-dpng','-r300')

%% Plot distance, silhouette and network templates
% cmap = distinguishable_colors(bestk);
cmap = CW.cMap; % hand-fixed networks
D_bestk = cell(bestk,1);
for inet = 1:bestk
    D_bestk{inet} = pdist2(template_X(assn_bestk==inet,:),C_bestk(inet,:),'correlation');
end

SI= silhouette(template_X,assn_bestk,'correlation');
mSI = mean(SI)

figure('position',[100 100 250 250]);
histogram(SI,-1:0.05:1,'EdgeColor','none','normalization','PDF');
xline(mSI,'Color','r');
xlabel('SI');
ylabel('PDF');
set(gca,'FontSize',12);
print(['./Figure/eLABE_template_bestk_SI_',method,'Zdim',num2str(zdim)],'-dpng','-r300');

figure('position',[100 100 800 400]);
xlim([0,bestk+1]);
ylim([-1,1]);
for inet = 1:bestk
    v = Violin(SI(assn_bestk==inet),inet,'ViolinColor',cmap(inet,:));
    v.ScatterPlot.SizeData = 10;
end
ylabel('SI');
xlabel('Networks');
xticks(1:bestk);
xticklabels(CW.Nets);
xtickangle(35);
set(gca,'FontSize',12);
print(['./Figure/eLABE_template_bestk_SI_violinplot_',method,'Zdim',num2str(zdim)],'-dpng','-r300');

N = ceil(sqrt(bestk));
figure('position',[100 100 800 800]);clear ax
t = tiledlayout(N,N,'TileSpacing','tight');
for inet = 1:bestk
    ax(inet) = nexttile;
    histogram(1-D_bestk{inet},'normalization','PDF','EdgeColor',[0.5 0.5 0.5],'FaceColor',cmap(inet,:),'FaceAlpha',1);
    set(gca,'FontSize',10);
    yticks([]);xticks([]);
    xlim([-1,1])
end
linkaxes(ax,'xy');
ylabel(t,'PDF','FontSize',12);
xlabel(t,'Correlation to centroid','FontSize',12);
print('./Figure/Correlation_to_bestk_centroids.png','-dpng','-r300');

% Across all networks
figure('position',[100 100 200 200]);
histogram(1-cat(1,D_bestk{:}),'normalization','PDF','EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',1);
ylabel('PDF');
xlabel('r');
xlim([-1,1])
set(gca,'FontSize',12);
print('./Figure/Correlation_to_bestk_centroids.png','-dpng','-r300');

%% Plot each template (centroid) and save them
cmap = ROY_BIG_BL(100);
colorrange = [-0.2,0.2];
for inet = 1:bestk
    vals = C_bestk(inet,:)';
    f = figure('position',[100 100 385 275]);
    ax1 = subplot(2,1,1);
    set(ax1,'Position',[0,0.5,0.85,0.5]);
    plot_parcels_by_values(vals,Anat,'lat',Parcels,colorrange,cmap)
    ax2 = subplot(2,1,2);
    set(ax2,'Position',[0,0.05,0.85,0.5]);
    plot_parcels_by_values(vals,Anat,'med',Parcels,colorrange,cmap)
    
    h = axes(f,'visible','off'); % attach colorbar to h
    c = colorbar(h,'Position',[0.88 0.1680 0.022 0.7],'XTick',[0,1],'XTicklabel',colorrange,'FontSize',12);
    colormap(c,cmap);
    exportgraphics(gcf,['./Figure/NetworkTemplate(bestk)_',num2str(inet),'.png'],'Resolution',300)
    close all
end

eval(['C_bestk',num2str(bestk),'=C_bestk;'])
eval(['assn_bestk',num2str(bestk),'=assn_bestk;'])
clear C_bestk assn_bestk
save([centroid_filename,'_centroids.mat'],'C_bestk*','assn_bestk*','-append');

