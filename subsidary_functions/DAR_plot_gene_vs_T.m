function [ output_args ] = DAR_plot_gene_vs_T( gene_names, thresh )
%	this function plots the Bjarnason data and REMAGUS data against the
%	estimated time. It colors the datapoints with Theta > thresh red.
%	Bjarnason data is green

times = {'am8','midday','pm4','pm8','midnight','am4'};
%load bjarnason data
load('frmaOM'); % Bjarnason data

% get Bjarnason times and Thetas
load('bjarnason_results.mat');
bjarnason_times=bjarnason.T;
bjarnason_Thetas=bjarnason.Thetas;

% get remagus data
load('esetfrmaBC.mat');

%get remagus results
load('remegus_results.mat');

load('Probes_string');

% get gene ids
gene_no = [];
index=[];
for i=1:length(Probes_string)
    for j=1:length(gene_names)
        str_probe = Probes_string{i}(1:length(gene_names{j}));
        str_gene = gene_names{j}(1:length(gene_names{j}));
        if strcmp(str_probe,str_gene)
             index = [index i];
             gene_no = [gene_no j];
        end
    end
end

for idx = 1: length(index)
    gene_id = index(idx);
    % put Bjarnason datainto array
    for tms = 1:6;
          for indivs=1:10;
            bjarnason_OMstr.(times{tms})(indivs) = frmaOM(gene_id,6*(indivs-1)+tms);
          end
    end
    bgenetimelist=zeros(2,60);
    % get bjarnason
    count = 1;
    for indiv=1:10
        for tms = 1:6
            genevect=bjarnason_OMstr.(times{tms})(indiv);
            T = bjarnason.T(tms,indiv);
            bgenetimelist(1,count) = genevect;
            bgenetimelist(2,count) = mod(T+6,24);
            count = count + 1;
        end
    end
    % get remagus
    rgenetimelist=zeros(2,226);
    count = 1;
    for i=1:226
        T = remegus.T(i);
        genevect=esetfrmaBC(gene_id,i);
        rgenetimelist(1,count) = genevect;
        rgenetimelist(2,count) = mod(T+6,24);
        rgenetimelist(3,count) = remegus.Theta(i);
        count = count + 1;
    end
    figure;
    scatter(rgenetimelist(2,:),rgenetimelist(1,:),'g');
    hold on
    scatter(bgenetimelist(2,:),bgenetimelist(1,:),'b');
    pcf=find(rgenetimelist(3,:)>thresh);
    scatter(rgenetimelist(2,pcf),rgenetimelist(1,pcf),'r')
    hold off
    end
end

