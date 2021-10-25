function Th_out = leave_1_out_Bjarn(excluded_probes,logthresh)

% This runs through individuals and removes the data for that individual
% from the trainingdata. It then uses the remining training data to
% cobstruct the probability model pm and then uses om
% calculate the times and Thetas for that individual whose data was excluded. 

% Out put is structure Th with fields
% D_Ts times calculated as Denise's thesis
% D_Thetas Thetas calculated as Denise's thesis
% Likelis_all 3d array with dimensions individual, time, projection time
% thus Likelis_all(individual,:,projection time) is the liklihood curve


times = {'am8','midday','pm4','pm8','midnight','am4'};

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

% make data structure
clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},'Affymetrix U133 2.0');

Th_out.D_Ts = zeros(length(times),10); % Denise's times
Th_out.D_Thetas = zeros(length(times),10); % Denise's Thetas

%Likelis_all =zeros(10,193,6);
% cycle through individuals
for individual = 1:10;
    %cs = clockOMstrN
    inds_data = [];
    % remove data for individual
	for tms = 1:6;
        %test(tms,:) = clockOMstrN2.(times{tms})(:,individual);
        cs = clockOMstrN.(times{tms});
        cs(:,individual)=[];
        cstr.(times{tms})=cs;
        inds_data = [inds_data clockOMstrN.(times{tms})(:,individual)];
    end
    % calculate probability model using the remining data
    pm = probabilty_model_construct(cstr,Inx,[1:9],{'am8','midday','pm4','pm8','midnight','am4'},0:4:20,24,7.5);
    % ge tthe likelihood curve for the individual's data
    [Likelis,RawLikelis] = get_likelis(pm,inds_data,[],logthresh);
    if individual == 1
        Likelis_all =zeros(10,size(Likelis,1),6);
    end
    % calculate times Ts and Thetas
    Th = Theta_calculator(Likelis);
    [D_Thetas, D_Ts] = D_calc_thetas(Likelis,pm);
%     Th.Likelis = Likelis;
%     Th.RawLikelis = RawLikelis;

    % put them into a structure
    Th_out.D_Ts(:,individual) = D_Ts'; % Denise's times
    Th_out.D_Thetas(:,individual) = D_Thetas'; % Denise's Thetas
    Th_out.Ts(:,individual) = Th.Ts; % Denise's times
    Th_out.Thetas(:,individual) = Th.Thetas; % Denise's Thetas
    Likelis_all(individual,:,:) = Likelis;
end
% put the likelihood curve into the structure
Th_out.Likelis_all=Likelis_all;
return

    