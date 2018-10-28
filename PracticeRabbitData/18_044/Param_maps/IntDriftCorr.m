%% Intensity drift correction 
cd('/v/raid2/sjohnson/Data/2018_RabbitData/18_044/Param_maps/'); 
load('Post_T1Map.mat'); 
post = paramImg; 
load('Pre_T1Map.mat');
pre = paramImg; 


figure; 
histogram(pre); 
hold on; 
histogram(post); 


for s = 1:size(post); 
    figure; histogram(pre(:,s,:),40,'BinLimits',[500 2500]); 
    hold on; histogram(post(:,s,:),40,'BinLimits',[500 2500]); 
    title(['slice = ' num2str(s)]); 
end 
%%
% Define the set of images as all slices in Pre and Post volumes 
sl_pre = permute(pre,[1,3,2]); 
sl_post = permute(post,[1,3,2]);

mask = zeros(size(sl_post)); 
for s = 1:size(sl_pre,3); 
    figure; 
    imagesc(squeeze(sl_pre(:,:,s)),[0 1500]); axis equal; 
    mask(:,:,s) = roipoly(); 
    close(gcf); 
end

sl_preM = sl_pre.*mask; 
sl_postM  = sl_post.*mask; 


Vpd = cat(3,sl_pre,sl_post);
Vpd_f = cat(3,sl_preM,sl_postM); 
numsl = size(sl_pre,3); 
%select a subset of images from Vpd representative of the range of
%histograms. For my case, this would be difference the top, bottom, and
%middle slices from Pre and Post volumes 
j = [2,ceil(numsl/2),numsl-2,numsl+2,numsl+ceil(numsl/2),2*numsl-2]
Vj = Vpd(:,:,j); 
Vj(Vj == 0) = NaN; 
maxV = max(Vj(:));

Vj_f = Vpd_f(:,:,j); 
Vj_f(Vj_f == 0) = NaN; 
maxV = max(Vj_f(:));

%%
figure; 
%for s = size(sl_preM,3); 
    
    V = squeeze(sl_preM(:,:,5)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on; 
    V = squeeze(sl_preM(:,:,10)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on; 
    V = squeeze(sl_preM(:,:,15)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on;
    V = squeeze(sl_preM(:,:,19)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on;
    
figure; 
%for s = size(sl_preM,3); 
    
    V = squeeze(sl_postM(:,:,5)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on; 
    V = squeeze(sl_postM(:,:,10)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on; 
    V = squeeze(sl_postM(:,:,15)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on;
    V = squeeze(sl_postM(:,:,19)) ;
    V = V(find(V >0)); 
    histogram(V,'Normalization','Probability'); hold on;
%end 
%%

% Load the pre-Map data, which has been resampled 
load('Pre_T1Map.mat'); 

%%

%define pc1 and pc2 and calculate the percentiles
pc1 =1; 
pc2 = 99;

m_1 = min(Vj(:))
m_2 = max(Vj(:))

%define [s1 s2] - which are the new mappings for pc1 and pc2
s1 = 1; %do not want to map to 0
s2 =4000; 
%%
for i = 4:6 %size(Vj,3); 
    Vi = Vj(:,:,i);
    Vi = reshape(Vi,[size(Vi,1)*size(Vi,2)],[]);
    
    % Make integer histogram and of the values in Vi
    intN = floor(m_1):ceil(m_2);
    N = histcounts(Vi,intN);
    %Mode_1 is the intensity value with the highest N counts between 800
    %and 1200
    N_sub = N(800-intN(1):1200-intN(1)); 
    intN_sub = intN(800-intN(1):1200-intN(1)); 
    [mN ind] = max(N_sub);
    m1 = intN_sub(ind);
    
    %calculate p1 and p2 from the perctentiles pc1 and pc2
    p1 = prctile(Vi,pc1);
    p2 = prctile(Vi,pc2);  
    
    figure; plot(intN(1:end-1),N);
    hold on; plot([p1,p2,m1],[0,0,0],'*','MarkerSize',8);
    
    %Map p1 and p2 onto s1 and 2 linearly 
    Vi_n = (Vi-p1).*((s2-s1)/(p2-p1))+s1; 
    
    %new histogram
    N_n = histcounts(Vi_n,intN); 
    hold on; plot(intN(1:end-1),N_n); 
    
    %calculate the new mode_1'
   % N_sub = N_n(s1:1500); 
   % intN_sub = intN(s1:1500); 
    [mN_n ind] = max(N_n);
    m1_prime = intN(ind);
    hold on; plot(m1_prime,0,'*','MarkerSize',8); 
    
    Mi(i) = m1; 
    P1i(i) = p1; 
    P2i(i) = p2; 
    Mi_pr(i) = m1_prime; 
end 

M_prmin = min(Mi_pr); 
M_prmax = max(Mi_pr); 
l = min(Mi-P1i); 
L = max(Mi-P1i); 
r = min(P2i-Mi); 
R = max(P2i-Mi);

% Condition2 term 
maxTerm = max([L/l,R/r]); 
cond = (L+R)*maxTerm; 

if (M_prmin-s1 < L) 
   % warning('Cond 1 for one-to-one mapping not met.');
    display(['[M_prmin - s1] ' num2str(M_prmin - s1) ' < [max(M - P1)] ' num2str(L) '.']); 
end 
if (s2-M_prmax < R)
   % warning('Cond 1 for one-to-one mapping not met.'); 
    display(['[M_prmax - s2] ' num2str(M_prmin - s2) ' < [max(P2 - M)] ' num2str(R) '.']);
end
if (s2-s1) < cond; 
   % warning('Cond 2 for one-to-one mapping not met.'); 
    display(['[(s2-s1)] ' num2str([s2-s1]) ' < [(L+R)*max(L/l,R/r)] ' num2str(cond) '.']); 
end

% check for one-to-one mapping


%%








load('/v/raid2/sjohnson/Data/2018_RabbitData/18_044/Seg3D_Masks/Day3_017_mrte_Muscle_CE_Final_MRTI_erode3.mat'); 
muscleVol = double(scirunnrrd.data); 

load('/v/raid2/sjohnson/Data/2018_RabbitData/18_044/TruncationInfo_forParameterMaps.mat'); 
truncMask = muscleVol(trR(1):trR(2),trC(1):trC(2),trS(1):trS(2)); 
postM = post.*truncMask; 
preM = pre.*truncMask; 

figure; 
%histogram(postM,'BinLimits',[1 2500]); 
histogram(postM,'BinLimits',[1 2500], 'Normalization','pdf'); 
hold on; 
histogram(preM,'BinLimits',[1 2500], 'Normalization','pdf'); 

%%
load('Post_T2Map.mat'); 
post = paramImg; 
load('Pre_T2Map.mat');
pre = paramImg; 
postM = post.*truncMask; 
preM = pre.*truncMask; 
figure; 
%histogram(postM,'BinLimits',[1 2500]); 
histogram(postM,'BinLimits',[1 100], 'Normalization','pdf'); 
hold on; 
histogram(preM,'BinLimits',[1 100], 'Normalization','pdf'); 

%% Make masks of hamstring muscle

mask = zeros(size(pre)); 
MRDataFig([],'pre',[0 100]); 
for s = 1:size(pre,2); 
    figure; 
    imagesc(squeeze(pre(:,s,:)),[0 100]); axis equal; 
    mask(:,s,:) = roipoly(); 
    close(gcf); 
end

preMham = pre.*mask; 
preMham = preMham(find(preMham > 0)); 
figure; histogram(preMham,'Normalization','pdf'); 
figure; histogram(preMham,'BinLimits',[600 1600]); 

postMham = post.*mask; 
postMham = postMham(find(postMham > 0)); 
hold on; histogram(postMham,'Normalization','pdf'); 
hold on; histogram(postMham,'BinLimits',[600 1600]);

