
function Th_out = leave_1_tissue_out(logthresh)
%
% This calculates the likelihood curves and Thetas for each of 8 tissues in the
% Zhang et al data using the probability model calculated using the other
% tissues.

load('zhang11_12organs_data')
% a 11x288 double array, 11 tissues

organnames={'Adrenal','Aorta','BrownFat','Heart','Kidney','Liver','Lung','SkelMus'};
times = {'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty','twentyone','twentytwo','twentythree','twentyfour'};

%clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');
for tms = 1:24; % times
    counter=0;
    for j=[1 2 4 6 8 9 10 11]; % 8 tissues
        counter = counter+1; % allows for possibility that some js are missing
        clockstr.(times{tms})(:,counter) = zhang11_12organs(:,24*(j-1)+tms);
        % clock11hogenstr.one is a 11x8 array, 11 probe expressions times 8 timepoints
    end
end
% gene vectors are normalised
for tms=1:24;
    for tis=1:8;
        clockstrN.(times{tms})(:,tis) = (clockstr.(times{tms})(:,tis)- mean(clockstr.(times{tms})(:,tis)))/std(clockstr.(times{tms})(:,tis));
    end
end
data_times=0:2:46;
instances = [1:8];
num_instances = length(instances);
pm=probabilty_model_construct(clockstrN,[],instances,times,data_times,48,3.75);
%%
Th_out.D_Ts = zeros(length(times),num_instances); % Denise's times
Th_out.D_Thetas = zeros(length(times),num_instances); % Denise's Thetas

Likelis_all =zeros(num_instances,length(pm.finetimes),length(data_times));
% cycle through individuals
for tissue = 1:8;
    %cs = clockOMstrN
    inds_data = [];
    % remove data for individual
	for tms = 1:24;
        %test(tms,:) = clockOMstrN2.(times{tms})(:,individual);
        % make clock str with the data for the tissue excluded
        cs = clockstrN.(times{tms});
        cs(:,tissue)=[];
        cstr.(times{tms})=cs;
        % put together data for the excluded tissue
        inds_data = [inds_data clockstrN.(times{tms})(:,tissue)];
    end
    % calculate probability model using the remining data
    pm=probabilty_model_construct(cstr,[],[1:8],times,data_times,48,3.75);
    % ge tthe likelihood curve for the individual's data
    %logthresh = -5;
    [Likelis,RawLikelis] = get_likelis(pm,inds_data,[],logthresh);
    % calculate times Ts and Thetas
    Th = Theta_calculator(Likelis);
    [D_Thetas, D_Ts] = D_calc_thetas(Likelis,pm);
%     Th.Likelis = Likelis;
%     Th.RawLikelis = RawLikelis;

    % put them into a structure
    Th_out.D_Ts(:,tissue) = D_Ts'; % Denise's times
    Th_out.D_Thetas(:,tissue) = D_Thetas'; % Denise's Thetas
    Th_out.Ts(:,tissue) = Th.Ts; % more accurate times
    Th_out.Thetas(:,tissue) = Th.Thetas; % more accurate  Thetas
    Likelis_all(tissue,:,:) = Likelis;
end
% put the likelihood curve into the structure
Th_out.Likelis_all=Likelis_all;