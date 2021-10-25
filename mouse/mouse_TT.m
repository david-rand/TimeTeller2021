function Th = mouse_TT( data_to_analyse,excluded_probes,periodicity,logthresh )
% This function takes the mouse dataset named data_str and calculates all the TT
% information from it.
% Use only for mouse data.
% Input: data_to_analyse: name of the data file e.g. 'esetfrmaBC' the REMAGUS data
% excluded_probes: any probes to be excluded from the rhythmic list
% periodicity use 48 or 24. The latter combines the training data at CT(t) with that at CT(t+24).
% The former analyses the data on a 48 hour basis to check the extent to
% which analysis of CT(t) and CT(t+24) agree.
% logthresh: sets the minimum level of the likelihood curves
%
%
% Outputs:
% D_Ts times using Denise's calculation
% D_thetas: Thetas using
% Ts more accurate times
% Thetas: more accurate Thetas
% Ratios: ratio between maximum of likelihood and the truncation level
% C_compare: the ratio of the rescaled truncation level (i.e. the truncation
% level of the likelihood ratio curve) and the minimum of C.

% Standard uses
% To apply TT to Bjarnason data
% [D_Ts, D_thetas, Ts, Thetas, Ratios, C_compare, Likelis, Log_max_likli, tmes] = human_TT(esetfrmaBC);

% Output structure:
% Th.D_Ts = D_Ts; % Denise's times
% Th.D_thetas = D_thetas; % Denise's Thetas
% Th.Ts = Ts; % more accurate times
% Th.Thetas = Thetas;% more accurate Thetas
% Th.Ratios = Ratios; % ratio of Log_max_likli to the logthreshold
% Th.C_compare = C_compare;
% Th.Likelis = Likelis; % the array of likelihoods
% Th.RawLikelis = RawLikelis; % the array of loglikelihoods before averaging
% Th.Log_max_likli = Log_max_likli; % the maximum value of the likelihood
% Th.tmes = tmes; the times correponding to finetimes, tmes(i) is time of
% ith point

% this is the code for leave one organ out validation

load('zhang11_12organs_data')
load('data_hogen')
Inx11=[29313 6747 25191 28920 7808 11341 14457 15237 29300 9574 26358];
% in case we need to exclude a gene as in Fang et al. Rev-Erb-Alpha is 11341
% which is probe 6, so to eclude this put excluded_probes = [6] in input
Inx=Inx11;
if ~isempty(excluded_probes)
    Inx(excluded_probes)=[];
    % Inx([6])=[];
end
% a 11x288 double array, 11 tissues and 
rhythmic_training_data = esetfrmaHogencore(Inx,:);
rhythmic_to_analyse = data_to_analyse(Inx,:);

Inx11names={"Arntl",
"Npas2",
"Per3",
"Dbp",
"Per2",
"Nr1d1",
"Nr1d2",
"Tef",
"Wee1",
"Ciart",
"Clock"};

organnames={'Adrenal','Aorta','BrownFat','Heart','Kidney','Liver','Lung','SkelMus'};
times = {'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty','twentyone','twentytwo','twentythree','twentyfour'};
times24 = {'thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty','twentyone','twentytwo','twentythree','twentyfour'};

%clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');
for tms = 1:24; % times
    counter=0;
    for j=[1 2 4 6 8 9 10 11]; % 8 tissues
        counter = counter+1; % allows for possibility that some js are missing
        clockstr.(times{tms})(:,counter) = rhythmic_training_data(:,24*(j-1)+tms);
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

% the above has periodicity 48 and uses all 24 times. For 24h periodicity
% we want to combine the training data at CT(t) with that at CT(t+24).
if periodicity == 24
    for tms=1:12;
        for tis=1:8;
            clockstrN.(times{tms})(:,tis+8)=clockstrN.(times24{tms})(:,tis);
        end
    end
    for i=13:24 clockstrN=rmfield(clockstrN,times{i});end
    data_times=0:2:22;
    pm=probabilty_model_construct(clockstrN,[],[1:16],times(1:12),data_times,24,15);
else
    pm=probabilty_model_construct(clockstrN,[],[1:8],times,data_times,48,15);
end
            

% pm=probabilty_model_construct(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');
%clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');

[Likelis,RawLikelis] = get_likelis(pm,rhythmic_to_analyse,[],logthresh);

[D_Thetas, D_Ts] = D_calc_thetas(Likelis,pm);

Th = Theta_calculator(Likelis);
Th.Likelis = Likelis;
Th.RawLikelis = RawLikelis;
Th.D_Ts = D_Ts; % Denise's times
Th.D_Thetas = D_Thetas; % Denise's Thetas
Th.logthresh = logthresh;
Th.Inx = Inx;

end

