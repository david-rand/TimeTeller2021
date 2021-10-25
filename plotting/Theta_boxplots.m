function Theta_boxplots(Th)

%Th_B = leave_1_out_Bjarn(excluded_probes,logthresh);

% box plot of individual's Thetas
groups=1:8;
groups=repmat(groups,8,1);
boxplot(Th.D_Thetas(:),groups(:));