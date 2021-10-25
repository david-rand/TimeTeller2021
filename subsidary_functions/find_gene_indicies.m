function indices = find_gene_indicies(gene_names)
% the names of the genes in the probe list are stired as 'genename_1053_at'
% e.g. 'RFC2_1053_at' and this function search for the genes in the list
% and returns their indices.

% Example: indices = find_gene_indicies({'PCNA'})

load('Probes_string');
gene_names_list = {};
indices=[];
% strcmp(Probes_string{i}(1:3),'HLF');
% gene_names= {'PCNA'}
% {'HLF' 'ARNTL2' 'NPAS2' 'TEF' 'RORC'};{'PCNA'}{'CKS2'}

gene_no = [];
for i=1:length(Probes_string)
    for j=1:length(gene_names)
        str_probe = Probes_string{i}(1:length(gene_names{j}));
        str_gene = gene_names{j}(1:length(gene_names{j}));
        if strcmp(str_probe,str_gene)
             indices = [indices i];
             gene_no = [gene_no j];
        end
    end
end
%%