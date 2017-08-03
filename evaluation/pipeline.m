
clear all;
close all;

% nohup matlab -nodesktop -nodisplay < er_community.m > mtestout.txt &

%% setup cvx
%cd cvx
%cvx_setup
%cd .. 

%% read in data, pre-process
chrnum = 20;
for cov = [10, 13, 17, 23, 26, 37]%2,4,6,8, 
disp(['Cov = ' num2str(cov)])
% adjacency matrix filled with +1/-1/0 
% final result is +1/+2
phase = importfile(['data/groundtruth/chr' num2str(chrnum) '_ground_truth.txt']);


% true ground truth file - consists of 1's and 2's
%ground_truth = zeros(1,length(phase));
%ground_truth(strcmp(phase,'0|1')) = 1;
%ground_truth(strcmp(phase,'1|0')) = 2;
%ground_truth(strcmp(phase,'0|2')) = 1;
%ground_truth(strcmp(phase,'2|0')) = 2;
% homozygous (multiple alternative alleles) - remove?
%ii1 = find(strcmp(phase,'1|2') == 1);
%ii2 = find(strcmp(phase,'2|1') == 1);
%ground_truth([ii1;ii2]) = [];   % 32778 
ground_truth = phase';
%length(ground_truth)




% adjacency matrix
%cd ../new_files/
data = load(['adjacent/chr' num2str(chrnum) 'adj_updated_cov' num2str(cov) '.csv'],'-ascii');
%cd ..


% index = []; ind = [];
% num = 1; num1 = 1;
% for i = 1:length(data)
%     if (abs(data(i,1)-data(i,2)) > 1500)
%         index(num) = i;
%         num = num+1;
%     end
% end
% data(index,:) = [];   % remove weird read that is far apart

%adj_num = length(data);
dat = spconvert(data); 

dat(end,end) = 0;   % only one with nonzero diagonal


% truncate to +1/-1/0 ?
dat(dat < 0) = -1;
dat(dat > 0) = 1;

% upper triangular all zeros, symmetrize
dat = dat+dat';

%size(dat)   % 32778*32778


% subsample the edges, needs to be symmetric
%cd groundtruth/index_files

%index = csvread('index_3.csv');
%new_dat = dat(index,index);
%data_sub = triu(dat,1);

%pp = find(data_sub~=0);
%pp1 = downsample(pp,5);   % subsampling ratio
%data_sub(setdiff(pp,pp1)) = 0;

%dat = data_sub+data_sub';
%size(dat)

% truncate ground truth file also
% ground_truth = ground_truth(index);
% length(ground_truth)

%cd ../../


%% construct moving & overlapping window for spectral analysis / SDP

W = 120;  % neighborhood size
%W = W - mod(W,4) + 2; 
%disp(['W = ' num2str(W)])
chunk_n = ceil((size(dat,2)-W)/(W/2));   % 655
mm = rem(size(dat,2)-W,W/2);  % number in last chunk  % 11
spec_res = zeros(chunk_n-1,W); 
err_rate = zeros(1,chunk_n); %err_rate_m = err_rate;
err_rate_sp = zeros(1,chunk_n);
spec_res_l = zeros(1,W);

tic
% multi-resolution scheme
% poolobj = parpool('local',2);
for i = 1:chunk_n-1   
    ww1 = (i-1)*W/2+1; ww2 = (i-1)*W/2+W;
    chunk = full(dat(ww1:ww2,ww1:ww2));   %W*W subgraph
    [S,C] = graphconncomp(sparse(chunk),'Directed','false','Weak','true');  % test if connected or not
    
    if S > 1  % has disconnected components

	tt = [];
	for q = 1:S
            tt = [tt find(C==q)];
	end        
	[~,I] = sort(tt);

    spec_res_l = [];
        for k = 1:S
            sub_chunk_ind = find(C==k)+(i-1)*W/2;
            sub_chunk = full(dat(sub_chunk_ind,sub_chunk_ind));
     
            [spec_tt,~] = spectral(sub_chunk,length(sub_chunk_ind),2,ground_truth(sub_chunk_ind));
            
            spec_res_l = [spec_res_l spec_tt];

        end

        spec_res(i,:) = spec_res_l(I);        
        
    else   % no disconnected componenets
        [spec_tt,~] = spectral(chunk,W,2,ground_truth(ww1:ww2));
        spec_res(i,:) = spec_tt;
    end
       
end


% last chunk
chunk = full(dat((chunk_n-1)*W/2+1:end,(chunk_n-1)*W/2+1:end));   % 110
[spec_res_fc,err_rate_sp(chunk_n-1)] = spectral(chunk,W+mm,2,ground_truth((chunk_n-1)*W/2+1:end));



%% put the result together - stitching

% final result vector
res_spec = zeros(1,size(dat,2));

% flipped version
spec_res_f = spec_res; %sdp_res_f = sdp_res; % MC_res_f = MC_res;
spec_res_f(spec_res==1) = 2;  spec_res_f(spec_res==2) = 1;

spec_res_uf = spec_res; % sdp_res_uf = sdp_res; % MC_res_uf = MC_res;


% overlap by W/2
for i = 1:chunk_n-2 
    
    % spectral method
    if sum(spec_res(i,W/2+1:W) == spec_res_uf(i+1,1:W/2)) < sum(spec_res(i,W/2+1:W) == spec_res_f(i+1,1:W/2))
        spec_res(i+1,:) = spec_res_f(i+1,:);   % flip 
        res_spec(i*W/2+1:(i+1)*W/2) = spec_res_f(i+1,1:W/2);
    else
        spec_res(i+1,:) = spec_res_uf(i+1,:); 
        res_spec(i*W/2+1:(i+1)*W/2) = spec_res_uf(i+1,1:W/2);
    end

end




% first chunk
res_spec(1:W/2) = spec_res_uf(1,1:W/2);


% last chunk
spec_res_fcf = spec_res_fc;% sdp_res_fcf = sdp_res_fc; % MC_res_fcf = MC_res_fc;
spec_res_fcf(spec_res_fc==1) = 2;  spec_res_fcf(spec_res_fc==2) = 1;


if sum(spec_res(chunk_n-1,W/2+1:W) == spec_res_fc(1:W/2)) < sum(spec_res(chunk_n-1,W/2+1:W) == spec_res_fcf(1:W/2))
    res_spec((chunk_n-1)*W/2+1:end) = spec_res_fcf;
else
    res_spec((chunk_n-1)*W/2+1:end) = spec_res_fc;
end




% check accuracy against ground truth - misclassification error

% nominal error
%accu_spec = sum(res_spec == ground_truth)/length(ground_truth)   %86  94
%accu_sdp = sum(res_sdp == ground_truth)/length(ground_truth)   %96  98


% switch error - better measure
% the proportion of SNPs whose phase is wrongly inferred relative to the
% previous SNP

% spectral method
tmp = zeros(1,length(res_spec)-1);
temp = zeros(1,length(res_spec)-1);
for i=1:length(res_spec)-1
    if res_spec(i) == res_spec(i+1)
        tmp(i) = 0;
    else
        tmp(i) = 1;
    end

    if ground_truth(i) == ground_truth(i+1)
        temp(i) = 0;
    else 
        temp(i) = 1;
    end
end
switch_err_spec = sum((tmp ~= temp))/length(tmp);

disp(['spec switch err: ' num2str(switch_err_spec)])

%delete(poolobj);

%disp('stitching done')

%%local correction step - Maximum Likelihood


%% initialize
res_spec_l = zeros(1,length(res_spec));

% only use local neighborhood reads for correction

%poolobj = parpool('local',2);
for refinenum = 1:3 % for parallel, uncomment previous line and use parfor
for tt = 1:length(res_spec)

    pp = find(dat(tt,:)~=0);
    pp(abs(pp-tt)>((refinenum)*W/2)) = [];    % remove long-range reads

    
    info = dat(tt,pp);
    known_spec = res_spec(pp);
    
    tinfo = info+known_spec;        
    final_spec = length(find(tinfo==3))+length(find(tinfo==0))-length(find(tinfo==2))-length(find(tinfo==1));
    
    if final_spec > 0
        res_spec_l(tt) = 2;
    else
        if final_spec < 0
            res_spec_l(tt) = 1;
        else
            res_spec_l(tt) = res_spec(tt);
        end        
    end

end

% update the result variable    
res_spec = res_spec_l;

    

% check accuracy - nominal error
%accu_spec = sum(res_spec == ground_truth)/length(ground_truth)  % 99


% switch error - better measure

% swe calculation for spectral method
tmp = zeros(1,length(res_spec)-1);
temp = zeros(1,length(res_spec)-1);
for i=1:length(res_spec)-1
    if res_spec(i) == res_spec(i+1)
        tmp(i) = 0;
    else
        tmp(i) = 1;
    end

    if ground_truth(i) == ground_truth(i+1)
        temp(i) = 0;
    else 
        temp(i) = 1;
    end
end
switch_err_spec = sum((tmp ~= temp))/length(tmp);

disp(['Round ' num2str(refinenum) ' spec switch err: ' num2str(switch_err_spec)]);

%disp('first round done')

end
toc

% save result
csvwrite(['data/res/chr' num2str(chrnum) '_cov' num2str(cov)],[res_spec;res_spec]);

end

