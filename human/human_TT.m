function Th = human_TT( data_to_analyse,excluded_probes,logthresh )
% This function takes the dataset named data_str and calculates all the TT
% information from it.
% Use only for human data.
% Input: data_str name of the data file e.g. 'esetfrmaBC' the REMAGUS data
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

load('frmaOM')
load('Inx16_Clean')
%load('Inx16')
load('Probes_string')
% datastr=load(data_str)
% names = fieldnames(datastr);

Inx = Inx16_Clean;
%Inx = Inx16;
if ~isempty(excluded_probes)
    Inx(excluded_probes)=[];
    % Inx([15])=[];
end

data_times = [0 4 8 12 16 20];

clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'});
pm=probabilty_model_construct(clockOMstrN,Inx,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},data_times,24,7.5);

 
% pm=probabilty_model_construct(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');
%clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');

[Likelis,RawLikelis] = get_likelis(pm,data_to_analyse,Inx,logthresh);

[D_Thetas, D_Ts] = D_calc_thetas(Likelis,pm);

Th = Theta_calculator(Likelis);
Th.Likelis = Likelis;
Th.RawLikelis = RawLikelis;
Th.D_Ts = D_Ts; % Denise's times
Th.D_Thetas = D_Thetas; % Denise's Thetas
Th.logthresh = logthresh;

end

