function probability_model = probabilty_model_construct(clockOMstr,rhythmic_probeset,instances_list,timeset,data_times,pm_period,pm_finetime)

% REPORT: Checked against OrigOrig_TimeTeller_REMAGUS 9/12/19 to get same
% pm using the Bjarnason data. Renamed prob_mod_construct.
% added

% OUTPUTS
% The output is a probability model structure. The fields are
%	Sigma: [1×1 struct]
%	Mu: [1×1 struct]
%	finetimes: [1×193 double]
%	proj: [1×1 struct]
%	rhythmic_probeset: [1×15 double]
%	probset: [54675×1 string]
%	original_times: {'am8'  'midday'  'pm4'  'pm8'  'midnight'  'am4'}
%	data_type: 'Affymetrix U133 2.0' : removed 20/12/19
%
% pm.Mu N d-dimensional vectors forming a circle in q-dimensional space
% each corresponding to a time t_i. These give the means of the multivariate
% normals. Each entry correaponds to a particlar time and N is closen so
% that the time gaps are small.
%
% pm.Sigma are the covariance matrices corresponding to the pm.Mu. They
% are periodic.
%
% pm.Times is anp array giving the times of each entry in pm.Mu and
% pm.Sigma
%
% pm.proj pm.proj(timeset{i}) gives the matrix M which defines the map from
% q-dimensional data space into d dimensions.
%
% pm.probset = probset this is the full probset supplied
%
% pm.rhythmic_probset = rhythmic_probset this is the found rhythmic probset
% that is used. The dimension q is length(pm.rhythmic_probset)
%
%
% Given data v from a microarray or similar we use the above to produce a
% likleihood function for v.

%  probabilty_model_construct(clockstr,rhythmic_probeset,instances_list,timeset,data_times,pm_period,pm_finetime)
% e.g.
% pm=probabilty_model_construct(cstr,Inx,[1:9],{'am8','midday','pm4','pm8','midnight','am4'},data_times,24,15);
% pm=probabilty_model_construct(clockstrN,[],[1:8],times,data_times,48,15);


% INPUTS
% clockOMstr: Data is input via clockOMstr which parameterised by instance (individual
% or tissue), probe and time. Thus clockstr.time1 is a num_probes x num_instances
% data set giving the gene expression for the chosen rhythmic genes for
% each instance. Only the rhythmic genes are included.
% rhythmic_probeset contains the gene ids of the rhythmic probes
% any dud probes should have been removed (eg probe 15 in the Bjarnason
% data)
% instances_list: means rhythmic_probes names of individual or tissues in the
% datasets we consider. For the 10 individuals in the Bjarnason data it is
% 1:10 and for the eight tssues in the Zhang et al data it is 1:8.
% timeset: timeset contains names for each time where data was collected in temporal
% order
% data_times are the numeric values corresponding to the times in timeset,
% pm_finetime is the time in minutesbetween times in the probability model
% pm_period is the period in hours used for the periodic time interval. This
% is usually 24 but can be 48. The latter is used for the Zhang moose data
% and is udes to compare the TT results at CT(t) with those at CT(t + 24).
% 




% % get rythmic gene names
% if ~isempty(rhythmic_probeset)
%     probenames = all_probe_names(rhythmic_probeset);
%     % % simplify them
%     probenames  = strrep(probenames , 'g', '');
%     probenames  =strrep(probenames , '_', ' ');
%     probenames  = strrep(probenames , 'at', '');
%     probenames  = strrep(probenames , 's', '');
%     probenames  = strrep(probenames , 'x', '');
%     % 
%     num_probes = length(rhythmic_probeset);
% end
num_times = length(timeset);
num_instances = length(instances_list);
% 
% % ===================== check compatibility of inputs ===================
% % =======================================================================
% 
% if size(expression_data,2) ~= num_times*num_instances
%     error('incompatability in inputs');
% end
% 
% % ========== make datastructure ==========
% % =======================================================================
% % make datastructure 
% rhythmic_data = expression_data(rhythmic_probeset,:); 
% % normalise data
% nrhythmic_data=(rhythmic_data-mean(rhythmic_data))./std(rhythmic_data);
% 
% % structure and normalise
% clear datastr
% for time_index = 1:num_times
% 	o = 0;
%     for i = 1:num_probes
%      % for i =[ 1 2 3 4 5 6 7 8 9 10 13 14];
%           o = o+1;
%           for j = 1:num_instances
%             datastr.(timeset{time_index})(o,j) = nrhythmic_data(i,num_times*(j-1)+time_index);
%           end
%     end
% end
% % datastr contains for each timepoint a matrix of 10 (num.
% % individs) numprobes-dimensional vectors that have been normalised to have
% % mean 0 and std 1.

datastr = clockOMstr;

% ========== do projections ==========
% =======================================================================
% for each time we calculate the projection into d dimensions
% first we SVD the data
clear A U V sig
d = 3;
for j = 1:num_times
    A.(timeset{j}) = datastr.(timeset{j})*(1/sqrt(5)); 
    [U.(timeset{j}),sig.(timeset{j}),V.(timeset{j})]=svd(A.(timeset{j})); 
end

% probability_model.XXX=A;
% return

% calculate the projections
% PROJECTION.(times{j}).(times{i})(1:3,:)' is projection into 3d
% The fields of PROJECTION are the times and, for example, PROJECTION.am8 has the
% same fields. Then e.g. PROJECTION.am8.midday is a 14×10 double obtined by
% using the am8 projection to project the midday data into R^3.
clear PROJECTION
for j = 1:num_times
	for i = 1:num_times
        PROJECTION.(timeset{j}).(timeset{i}) = U.(timeset{i})'*datastr.(timeset{j}); 
            % this is the data for time t_j projectred using the projection
            % from time t_i
    end
end
% probability_model.XXX=PROJECTION;
% return
% ========== construct means and covariances ==========
% =======================================================================
% GMModel = fitgmdist(X,k,Name,Value) returns a Gaussian mixture distribution model 
% with additional options specified by one or more Name,Value pair arguments.
% Then e.g. GMModel.am8.midday is a fit of a multivariate normal to
% PROJECTION.am8.midday in 3d
% PROJECTION.am8.midday uses the am8 projection to project the midday data into R^3 and fit a normal distribtion.
clear GMModel
for j = 1:num_times
    for i = 1:num_times
        GMModel.(timeset{j}).(timeset{i}) = fitgmdist(PROJECTION.(timeset{j}).(timeset{i})(1:d,:)',1);
    end
end


% ========== Now we interpolate the distributions to all times ==========
% =======================================================================
% pm_period = 24;
% data_times=[0, 4, 8, 12, 16, 20]';
% periodic_data_times = [0, 4, 8, 12, 16, 20, 24]';
% 
% 
% num_cycles = 1;
% time_spacing = 4;
% time_length = 24;
% int = 0.125/2; % interval between data points



signam = {'one','two','three','four','five','six','seven','eight','nine', 'ten'};

fine_times = 0:(pm_finetime/60)*(pm_period/48):pm_period; % all ellipsoids. This is a 15 min interval
data_time_step = data_times(2)-data_times(1);
data_periodic_times = [0:data_time_step:length(data_times)*data_time_step];

for projtime = 1:num_times % projecting time
    clear mus sigmas pp ppm mus_all Sigma_mat 
    for ti = 1:num_times
        mus(:,ti) = GMModel.(timeset{ti}).(timeset{projtime}).mu;
        % mus(:,j) is the mean of the t_j data points after projection by
        % the t_var projector
        sigmas(:,ti) =   GMModel.(timeset{ti}).(timeset{projtime}).Sigma(:);
        % as above but covariance matrices
    end
    mus(:,num_times+1) = GMModel.(timeset{1}).(timeset{projtime}).mu;
    % this makes  mus(:,7) = mus(:,7)
    sigmas(:,num_times+1) = GMModel.(timeset{1}).(timeset{projtime}).Sigma(:);
    for i = 1:d^2
        pp1.(signam{i}) = csape(data_periodic_times',[sigmas(i,:)'],'periodic');
        % can also use pp2.(signam{i}) =
        % spcsp(data_periodic_times',[sigmas(i,:)],w) with p a smoothing
        % parameter where csape correponds to p=1, cam also weight vectors
        % pp = pchip(x,y) returns a piecewise polynomial structure for use with ppval and the spline utility unmkpp.
        % p=pchip(x,y,xq) returns a vector of interpolated values p 
        % corresponding to the query points in xq. The values of p 
        % are determined by shape-preserving piecewise cubic interpolation of x and y.
    end
    % matches first and second derivatives at end points
    all_mus=mus(1,:)';for i=2:d all_mus=[all_mus,mus(i,:)'];end
    ppm = csape(data_periodic_times',all_mus','periodic');
    mus_all = ppval(ppm,fine_times); 
    % mus_all is 3 x N where N is the number of interpolating points: it
    % contains the interpolation of the means.
    % when d=3 all_mus = [mus(1,:)',mus(2,:)',mus(3,:)']' is a periodic array of
    % num_time + 1 d-dimensional vetcors

    for i = 1:d^2
        Sigma_mat(i,:) = ppval(pp1.(signam{i}),fine_times); 
        % this is a d^2 x N array each of the N entries a dxd matrix
    end

    % checking for non positive eigenvalues (used rarely, but necessary)
    countq = 0;
    for q = 1:length(mus_all);
        if find(eig(reshape(Sigma_mat(:,q),d,d)) <=0.0001)
            A=reshape(Sigma_mat(:,q),3,3);
            B = (A + A')/2; 
            [U5,Sigma5] = eig(B); %
            Ahat5 = U5*max(Sigma5,0)*U5'; %
            Ahat5=Ahat5+0.001*eye(3);
            Sigma_mat(:,q) =Ahat5(:);
            countq = countq+1;
        end
        % put the mus and Sigmas into a structure indexed by the projecting
        % time and the index q in the interpolating timeseries
        Sigma_struct.(timeset{projtime})(:,q) = Sigma_mat(:,q);
        Mu_struct.(timeset{projtime})(:,q) = mus_all(:,q);
    end
    % Sigma_struct and Mu_struct give the MVN distribution for all times.
end

probability_model.Sigma = Sigma_struct;
probability_model.Mu = Mu_struct;
probability_model.finetimes = fine_times;
for i = 1:length(timeset)
     probability_model.proj.(timeset{i}) = U.(timeset{i})';
     probability_model.U.(timeset{i}) = U.(timeset{i});
end
if ~isempty(rhythmic_probeset)
    probability_model.rhythmic_probeset = rhythmic_probeset;
    % probability_model.probset = all_probe_names; 
end
probability_model.original_times = timeset;
probability_model.datastr = datastr;
    
return

