function  datastr = make_clock_OM_str(exp_data,rhythmic_probeset,all_probe_names,instances_list,timeset,data_type)

% This takes in the Bjarnason data or similar data and puts it in arrays
% indexed by time. Each of these arrays is (data dimension)x(no. of
% instaces). For the Bjarnson data the instances are the ten individuals.

% The data structire datastr has fields given by timeset. This is usually
% {'am8','midday','pm4','pm8','midnight','am4'}.

expression_data = exp_data;

%show('Have you removed dud probes?')

% get rythmic gene names
probenames = all_probe_names(rhythmic_probeset);
% simplify them
probenames  = strrep(probenames , 'g', '');
probenames  =strrep(probenames , '_', ' ');
probenames  = strrep(probenames , 'at', '');
probenames  = strrep(probenames , 's', '');
probenames  = strrep(probenames , 'x', '');

num_probes = length(probenames);
num_times = length(timeset);
num_instances = length(instances_list);

% ===================== check compatibility of inputs ===================
% =======================================================================

if size(expression_data,2) ~= num_times*num_instances
    error('incompatability in inputs');
end

% ========== make datastructure ==========
% =======================================================================
% make datastructure 
rhythmic_data = expression_data(rhythmic_probeset,:); 
% normalise data
nrhythmic_data=(rhythmic_data-mean(rhythmic_data))./std(rhythmic_data);

% structure and normalise
clear datastr
for time_index = 1:num_times
	o = 0;
    for i = 1:num_probes
     % for i =[ 1 2 3 4 5 6 7 8 9 10 13 14];
          o = o+1;
          for j = 1:num_instances
            datastr.(timeset{time_index})(o,j) = nrhythmic_data(i,num_times*(j-1)+time_index);
          end
    end
end
% datastr contains for each timepoint a matrix of 10 (num.
% individs) numprobes-dimensional vectors that have been normalised to have
% mean 0 and std 1.
return

% 
% This function turns the Bjarnason data into clockOMstrN
% where  clockOMstr.(times{tmes})(o,tms) is 
% estimated times, (iii) plots the likelihood functions, 
% (iv) plots a histogram of Thetas, (v) plots more likelihood functions: 
% set up verticals and then plot likelihoods.

% History
% 2016 Used for Denise Thesis
% 2018-9 Used for paper.
% Last check: May 2019
% Further commenting: May 2019

disp('Running make_clock_OM_str');

load('Inx16') % indexes of most synchronised rhyhtmic probes
%load('Probes_cell')
load('Probes_string')
load('frmaOM')

times = {'am8','midday','pm4','pm8','midnight','am4'};

cm=hsv(6);

% the full set of genes found by the synchronicity and rhythmicity analysis
% have the probe number in Inx16 
Inx=Inx16;
% with names given by
probenames = Probes_string(Inx16);
%just making the names easier to read in figures
probenames  = strrep(probenames , 'g', '');
probenames  =strrep(probenames , '_', ' ');
probenames  = strrep(probenames , 'at', '');
probenames  = strrep(probenames , 's', '');
probenames  = strrep(probenames , 'x', '');

% We use 16 probes here and the corresponding data from Bjarnason et al 
% is contained in bjarn16
load('Bjarn16');
clear clock_OM
clock_OM = bjarn16;


num_times = length(timeset);
num_instances = length(instances_list);

% %manually change
% clear clock_OM
% Inx=Inx16; 
% %% comment out for histogram, leave in for LOOV
Inx([15])=[]; %% here we delete the 15th probe based on the observation that probe intensites are inconsistent amonst independent experiements
clock_OM = frmaOM(Inx,:); % picking the 15 circadian timecourses (this can be saved and loaded for publishing purposes)
num_probes = length(Inx);
%% structure and normalise
clear clockOMstr
for tmes = 1:6; % number of timepoints
    g=0; %genes
    for i = 1:length(Inx); %list of probes
        g = g+1;
        if i ~= g disp('i not equal g');end
        for tms=1:10; % individuals
            clockOMstr.(times{tmes})(g,tms) = clock_OM(i,6*(tms-1)+tmes);
        end
    end
end

% for tmes = 1:6; % number of timepoints
%     for tms=1:10; % individuals
%     clockOMstr.(times{tmes})(:,tms) = clock_OM(:,6*(tms-1)+tmes);
%     end
% end

% make datastructure 
rhythmic_data = clock_OM; 
% normalise data
nrhythmic_data=(rhythmic_data-mean(rhythmic_data))./std(rhythmic_data);

for time_index = 1:num_times
	o = 0;
    for i = 1:num_probes
     % for i =[ 1 2 3 4 5 6 7 8 9 10 13 14];
          o = o+1;
          for j = 1:num_instances
            datastr.(timeset{time_index})(o,j) = nrhythmic_data(i,num_times*(j-1)+time_index);
          end
    end
end
clockOMstrN=datastr;
% clockOMstr
% clear clockOMstrN
% for tmes=1:6;
% 	for indivs=1:10;
%         v = clockOMstr.(times{tmes})(:,indivs);
%         v = (v - mean(v))/std(v);
%         clockOMstrN.(times{tmes})(:,indivs) = v;
% 	end
% end
return

