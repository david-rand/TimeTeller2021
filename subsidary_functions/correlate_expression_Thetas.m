function [r,p]=correlate_expression_Thetas(gene_index, Thetas, corr_type)

% calculates correlation between Thetas and the gene with the given index.
% corr_type can be 'Pearson', 'Kendall' or 'Spearman'

% note index of PCNA is 10651
% 13618 is index of CKS2

load('esetfrmaBC')
expressions = esetfrmaBC(gene_index,:);
figure
scatter(expressions,Thetas)
[r,p]=corr(expressions',Thetas','Type',corr_type);
return