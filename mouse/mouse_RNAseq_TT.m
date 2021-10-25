function Th = mouse_RNAseq_TT( data_to_analyse,excluded_probes,periodicity,logthresh )
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

load('rnaseqzhang')
load('newTop15Genes')
% a 11x288 double array, 11 tissues

rhythmic_to_analyse = data_to_analyse(newTop15Genes,:);

clk_rnaseqzhang=rnaseqzhang(newTop15Genes,:);

organnames={'Adrenal','Aorta','BrownFat','Heart','Kidney','Liver','Lung','SkelMus','WhiteFat'};
times = {'one','two','three','four','five','six','seven','eight'};

%clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');
for tms = 1:8; % times
    counter=0;
    for j=2:9; % 8 tissues
        counter = counter+1; % allows for possibility that some js are missing
        clockstr.(times{tms})(:,counter) = clk_rnaseqzhang(:,8*(j-1)+tms);% 8 is no timepoints
        % clock11hogenstr.one is a 11x8 array, 11 probe expressions times 8 timepoints
    end
end
% gene vectors are normalised
for tms=1:8;
    for tis=1:8;
        clockstrN.(times{tms})(:,tis) = (clockstr.(times{tms})(:,tis)- mean(clockstr.(times{tms})(:,tis)))/std(clockstr.(times{tms})(:,tis));
    end
end
data_times=0:6:42;
instances = [1:8];
num_instances = 8;

times = {'one','two','three','four','five','six','seven','eight'};
times24 = {'five','six','seven','eight'};
% the above has periodicity 48 and uses all 24 times. For 24h periodicity
% we want to combine the training data at CT(t) with that at CT(t+24).
if periodicity == 24
    for tms=1:4;
        for tis=1:8;
            clockstrN.(times{tms})(:,tis+8)=clockstrN.(times24{tms})(:,tis);
        end
    end
    for i=5:8 clockstrN=rmfield(clockstrN,times{i});end
    data_times=0:6:18;
    pm=probabilty_model_construct(clockstrN,[],[1:16],times(1:4),data_times,24,15);
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
Th.Inx = newTop15Genes;

end

