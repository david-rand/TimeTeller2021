function Th = Theta_calculator(lik)

% last amended 3/12/19 by DAR, 2/1/2022 by DAR
% this is a more accurate way of finding times and Thetas. It uses splines
% to calculate the crossings of C and the likelihood ratio curves more
% accurately.
% lik is an array of likelihood curves each of which is a column vector with equal
% time spacing

% Output structure:
% Th.Ts = Ts; % more accurate times
% Th.Thetas = Thetas;% more accurate Thetas
% Th.Ratios = Ratios; % ratio of Log_max_likli to the logthreshold
% Th.C_compare = C_compare;
% Th.Log_max_likli = Log_max_likli; % the maximum value of the likelihood
% Th.tmes = tmes; the times correponding to finetimes, tmes(i) is time of
% ith point
% Th.second_Ts gives the time of the second peak if there is one i.e. the
% highest peak that is at least 6h away from the main peak. If curve is
% flat gives NaN, if only one peak gives -1.
% Th.second_Ts_ratio gives the ratio of the likelihood of the second peak
% to that of the main one

% Parameters for C
epsilon=0.4;
eta=0.35;

% time points scaled to be between 0 and 1
x = (1:size(lik,1))/size(lik,1);x=x';
tmes =24*x;
tmes = mod(tmes+8,24);


Ts = []; % will contain times on (0,1)
Tclock = [];% real time assuming Ts = 0 corresponds to 8am
Ratios = []; % will contain ratio between flat level and max of likleli
C_compare = [];% will contain ratio of Ratio to min of C: less than 1 implies no intersection
Thetas = zeros(1,size(lik,2));% will contain Thetas
% ddots = [];
goodindices = [];
badindices = [];
indices= [];
Log_max_likli = [];
count_flats = 0;
flats = [];
sT=[];
second_ratio = [];
num_intersections_list = [];

% loop through the likelihoods
for num_likeli = 1:size(lik,2)
    y=lik(:,num_likeli);   
    I = find(y>min(y))'; % indicies outside flat regions
    if ~isempty(I)
        J=I;J(1)=[];K=J-I(1:end-1);K=[2 K];sts=find(K>1);eds = sts-1;eds=eds(2:end);eds = [eds length(I)];
        num_nonflats = length(sts);

        %make spline
        pp=csape(x,y,'periodic');
        err = max(abs(y-ppval(pp,x)));
        criticals =[];
        for i=1:num_nonflats %loops over nonflat regions and finds all critical points in them
            if x(I(sts(i)))~= x(I(eds(i)))
                criticals = [criticals fnzeros(fnder(pp),[x(I(sts(i))),x(I(eds(i)))])];
%             else
%                 criticals = [criticals fnzeros(fnder(pp),[x(I(sts(i))-1),x(I(eds(i))+1)])];
            end
        end

        if criticals(1,:) ~= criticals(2,:)
            error('criticals do not match')
        end
        criticals = criticals(1,:); % list of critical points

        z=ppval(pp,criticals);
        [m ind]=max(z(1,:));
        T=criticals(ind); % find time T
        Ts = [Ts T];
        Tclock = [Tclock mod(24*T+8,24)];
        r = min(y)/m;
        Ratios = [Ratios, r];
        Log_max_likli = [Log_max_likli log(m)];

        if isempty(T)
            error('T empty')
        end
        if T > 1 T = T - 1; end
        lr = y/m;
        C=1+epsilon + cos(2*pi*(x-T));
        C=eta*C;
        C_compare = [C_compare r/min(C)];

        pp2=spline(x,lr-C); % fit spline
        zros = fnzeros(pp2); % find zeros
        zros = zros(1,:);
        num_intersections = length(zros);
        zros = [0 zros 1]; % add endpoints to zeros
        len = 0;
        for k = 1:length(zros)-1
            u =  zros(k) + (zros(k+1)-zros(k))/2;
                % check midpoint of each region to see if lr-C is positive or
                % negative
                if ppval(pp2,u)>0
                    len = len + zros(k+1)-zros(k); %add the ones where lr-C is positive
                end
        end
        Thetas(num_likeli) = len; % this gives Theta
        indices=[indices num_likeli];
        % secondpeak stuff
        if 1%mod(24*T+8,24)<10%secondpeak
            %sec_criticals=[];
            %plot(x,y)
            Tl=T-0.25;
            Tr = T+0.25;
            if Tl < 0
                sec_criticals = criticals(find(and(criticals>T+0.25,criticals<T+0.75)));
            else if Tr >1
                sec_criticals = criticals(find(and(criticals>T-0.75,criticals<T-0.25)));
                else
                    sec_criticals = criticals(find(or(criticals<Tl,criticals>Tr)));
                end
            end
            z2=ppval(pp,sec_criticals);
            if ~isempty(z2)
                [m2 ind2]=max(z2(1,:));
                sT=[sT sec_criticals(ind2)];
                %plot(x,y)
                %disp(sec_criticals)
                %disp(sT)
                second_ratio=[second_ratio m2/m];
            else
                sT=[sT -1];
                second_ratio=[second_ratio -1];
            end
                
        end
    else % flat curve
        Thetas(num_likeli) = 1;
        num_intersections = 0;
        Ratios = [Ratios, 0];
        indices=[indices num_likeli];
        Ts = [Ts -1];
        Tclock = [Tclock -1];
        Log_max_likli = [Log_max_likli -Inf];
        count_flats = count_flats + 1;
        flats = [flats num_likeli];
        sT=[sT NaN];
        second_ratio=[second_ratio NaN];
    end
    num_intersections_list = [num_intersections_list num_intersections];
end
Th.Ts = Ts; % more accurate times
Th.Thetas = Thetas;% more accurate Thetas
Th.Ratios = Ratios; % ratio of Log_max_likli to the logthreshold
Th.C_compare = C_compare;
Th.Log_max_likli = Log_max_likli; % the maximum value of the likelihood
Th.tmes = tmes;
Th.Tclock=Tclock;
Th.flats = flats;
Th.second_Ts = sT;
Th.second_Ts_ratio=second_ratio;
Th.num_intersections = num_intersections_list;

return