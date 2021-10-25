function [Thetas, D_Ts] = D_calc_thetas(M2,pm)
% 9/12/19
% This implements Denise's method of calculating Theta.
% M2 is an array of likelihoods represented by the columns and pm is the
% probability model which was used to calculate the likelihoods

finetimes = size(M2,1);
num_vects = size(M2,2);
times = pm.original_times;


% for p = 1:num_vects;
%     [L(p), a(p)] = max(M(:,p));
% end
 
for p = 1:num_vects;
    [L2(p), a2(p)] = max(M2(:,p));
end
 
% a = a/finetimes *24;
a2 = a2/finetimes *24;
 
% for r =1:6;
%     lengths(r) = length(pm.Mu.(times{r}));
% end
 
% for q=1:min(lengths);
%     %M(q) = mean(m3(:,q));
%     M2(q) =  geomean(m4(:,q));
% end
 
% parameters for C
eta=0.35; 
epsil=0.4;

halftimes = floor(finetimes/2); %96
 
clear dys_new prob prop Forward Backward LFN LBN LF LB
for yu = 1:num_vects;
%figure(); plot(M2(:,yu))
    T = find(M2(:,yu)== max(M2(:,yu)));
    T=T(1);
    D_Ts(yu) = T;
    if T>halftimes
        Forward = [T+1:finetimes,1:(T-halftimes-1)];
        Backward = [T-halftimes:T-1];
    else
         Forward = [T+1:(T+halftimes)];
        Backward = [halftimes+1+T:finetimes,1:T-1];
    end
    %% check
    %figure(1); hold on;  plot(Forward,M2(Forward,yu));
     %hold on; plot(Backward,M2(Backward,yu)) 
    for c =1:halftimes;
    LF(c) = 1/(epsil+1+cos(c/halftimes*pi)).*M2(Forward(c),yu);
    end
    for c =1:halftimes;
    LB(c) =1/(epsil+1+cos(((halftimes+1)-c)/halftimes*pi)).*M2(Backward(c),yu);
    end

    % LFN = M2(T,yu)./LF;
    % LBN = M2(T,yu)./LB;

    LFN = LF./M2(T,yu);
    LBN = LB./M2(T,yu);

    %figure(2); plot(Backward,LBN,'*'); hold on; plot(Forward,LFN,'*')
    prop = find([LFN LBN]>eta);
    prob(yu) = length(prop)/finetimes; 
    %dys_new(yu) =log10(sum([LFN LBN]));
end
Thetas = prob;
return