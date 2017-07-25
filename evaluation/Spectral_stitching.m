function [res_spec] = Spectral_stitching(Input_file, Output_file)


% nohup matlab -nodesktop -nodisplay < er_community.m > mtestout.txt &

%% setup cvx
%cd cvx
%cvx_setup
%cd .. 

%% read in data, pre-process


% adjacency matrix filled with +1/-1/0 
% final result is +1/+2

data = load(Input_file,'-ascii');

dat = spconvert(data); 

dat(end,end) = 0;   % only one with nonzero diagonal


% truncate to +1/-1/0 ?
dat(dat < 0) = -1;
dat(dat > 0) = 1;

% upper triangular all zeros, symmetrize
dat = dat+dat';

% subsample the edges, needs to be symmetric


%% construct moving & overlapping window for spectral analysis / SDP

W = 120;  
chunk_n = ceil((size(dat,2)-W)/(W/2));   % 655
mm = rem(size(dat,2)-W,W/2);  % number in last chunk  % 11
spec_res = zeros(chunk_n-1,W); 
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
     
            [spec_tt,~] = spectral(sub_chunk,length(sub_chunk_ind),2,-1);
            
            spec_res_l = [spec_res_l spec_tt];

        end

        spec_res(i,:) = spec_res_l(I);        
        
    else   % no disconnected componenets
        [spec_tt,~] = spectral(chunk,W,2,-1);
        spec_res(i,:) = spec_tt;
    end
       
end


% last chunk
chunk = full(dat((chunk_n-1)*W/2+1:end,(chunk_n-1)*W/2+1:end));   % 110
[spec_res_fc,err_rate_sp(chunk_n-1)] = spectral(chunk,W+mm,2,-1);



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

    
disp(['Round ' num2str(refinenum) '  refinement finished.']);

end
toc

% save result
csvwrite(Output_file, res_spec);

end

