function [Likelis,RawLikelis] = get_likelis(pm,expression_data,rhythmic_probeset,logthresh)
% Output
% pm is probability model to be used
% Inputs:
% expression_data is the expresseion data to be analkysed (e.g. frmaOM)
% rhythmic_probeset are the rythmic data indices (e.g. Inx16 with 15 removed)
% Example: Thetas = get_Thetas(pm,frmaOM,Inx)

% ========== make datastructure ==========
% =======================================================================
% make datastructure 
if ~isempty(rhythmic_probeset)
    rhythmic_data = expression_data(rhythmic_probeset,:); 
else
    rhythmic_data = expression_data;
end
% rhythmic_data = expression_data;
num_vects = size(rhythmic_data,2);
% normalise data
nrhythmic_data=(rhythmic_data-mean(rhythmic_data))./std(rhythmic_data);

% probability_model.proj.(timeset{i}) = U.(timeset{i})';
% PROJECTIONN.test.(times{i}) = UN.(times{i})'*BCN;
times = pm.original_times;
Utranspose = pm.proj;
Mu_struct = pm.Mu;
Sigma_struct = pm.Sigma;

% These are the projections of the data onto the first d principal
% components
for i=1:length(times)
	PROJ.(times{i}) = Utranspose.(times{i})*nrhythmic_data;
end

% this calcul;ates the likelihoods q indexes time
for r = 1:length(times);
    for q =1:length(pm.finetimes);
        m(r,q,:) = log(mvnpdf(PROJ.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3)));  
        %m2(r,q,:) = mvnpdf(PROJ.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3));
    end
end
m3=m;  
% m4 = m2;
% logthresh = -10;
explogthresh = exp(logthresh);

% Apply the threshold to each likeli and average ocross the different
% projection times
for p =1:num_vects;
    for q=1:length(pm.finetimes);
        for r = 1:length(times);
            m3(r,q,p) = max(m3(r,q,p),logthresh); %-5,-10
            %m4(r,q,p) = max(m4(r,q,p),explogthresh);
        end
    M(q,p) = sum(m3(:,q,p))/length(times);
    %M2(q,p) =  geomean(m4(:,q,p));
    M2 = exp(M);
    end
end

% max(max(abs(M2-exp(M))))

Likelis = M2;
RawLikelis = m;
% MaxLikelis=squeeze(max(m));