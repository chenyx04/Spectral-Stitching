

% identify cluster using spectral method 
function [idx,err_rate] = spectral(We,n,Q,truth)

    [V,D] = eig(We);
    eval = diag(D);
    [~,I] = sort(eval,'descend');
    V1 = V(:,I);
    evec = V1(:,1:Q-1);  
 
    % thresholding the top eigenvector 
    idxt = ones(1,n);
    idxt(evec >= 0) = 2;
    idx = idxt;
%     [~,idx] = compute(idxt,truth);
if truth~=-1
    % calculate switch error
    tmp = zeros(1,length(idx)-1);
    temp = zeros(1,length(idx)-1);
    for i=1:length(idx)-1
        if idx(i) == idx(i+1)
            tmp(i) = 0;
        else
            tmp(i) = 1;
        end
        
        if truth(i) == truth(i+1)
            temp(i) = 0;
        else 
            temp(i) = 1;
        end
    end
    err_rate = sum((tmp ~= temp));
else
    err_rate=0;
end
%     % run 5 times with different initialization for K-means
%     idx_t = zeros(5,n);  val = zeros(1,5);
%     for tt = 1:5
%         idx_t(tt,:) = kmeans(evec,Q).';
%         [val(tt),idx_t(tt,:)] = compute(idx_t(tt,:),truth);   % permutation of labels
%     end
%     
%     [~,pick] = min(val);
%     idx = idx_t(pick,:);

    
%     idx = kmeans(evec,Q,'Replicates',5).';
    
return




