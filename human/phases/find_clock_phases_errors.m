function find_clock_phases_errors(excluded_probes)

% Input excluded_probes, a list of rthymic probes to be excluded e.g. [11 15].
% For this routine usually no probes are exluded

% History.
% Written 20/12/19 using human_TimeTeller_phases from Denise
% Tested 20/12/19

% This function calculates the cosinor phase for each Bjarnason
% data point and then calculates T and the likelihood curve for them 
% using a leave one out approach. It then plots (i) the residual (x-axis) 
% against the Cosinor phase of the gene (y-axis), (ii) the likelihood
% curves and (ii) the actual time vs estimated time in a scatter plot

% there are a total of 60 datapoints = 6 timepoints x 10 individualks

% This uses all 16 probes

load('frmaOM') % is this data meddled with? or was the data before looking the way it was due to a shift in the data?
load('frmaOM17_individual_labels')
load('clock16OMJULY')
load('clock16OMnames')
load('Inx16_Clean')
%load('Inx16')
load('Probes_string')
% datastr=load(data_str)
% names = fieldnames(datastr);

times = {'am8','midday','pm4','pm8','midnight','am4'};

Inx = Inx16_Clean;
%Inx = Inx16;
if ~isempty(excluded_probes)
    Inx(excluded_probes)=[];
    % Inx([15])=[];
end

num_probes = length(Inx);

cm=hsv(6);
probenames = Probes_string(Inx); % takes all probstrings and restricts to only those used
probenames  = strrep(probenames , 'g', '');
probenames  =strrep(probenames , '_', ' ');
probenames  = strrep(probenames , 'at', '');
probenames  = strrep(probenames , 's', '');
probenames  = strrep(probenames , 'x', '');
% probenames lists the 16 genes/probes being used - the steps are just
% removing unwanted parts of the stored name strings from 

%manually change
clock_OM = frmaOM(Inx,:); %the vectors for the 60 datasets

bmal=clock_OM(7,:); %data for  Bmal1

Th_out = leave_1_out_Bjarn([],-10);

% use cosinor to get the phases of the probes for each individual
clear q a p allphases maxf minf antiphase
for q=1:num_probes; % runs over probes
    for a = 1:10; % runs over individuals
        [p(q,a) allphases(q,a),maxf(q,a),antiphase(q,a), minf(q,a)]= cosinor(linspace(0,1,7),clock_OM(q,[(a-1)*6+1:6*a, (a-1)*6+1]),2*pi,0.05);
        if allphases(q,a)>24
            allphases(q,a)=allphases(q,a)-24;
        end
    end
end

a2=Th_out.D_Ts;
%a2 = a2/193 *24;
a2 = mod(24*a2/193+8,24);
%actuals = repmat(0:4:20,10,1)';
actuals = repmat([8 12 16 20 0 4],10,1)';

errors = actuals-a2; %difference estimated time and actual time
% must allow for fact that 24h error is not error
for t = 1:60;
    if errors(t)>12
        errors(t)=errors(t)-24;
    end
    if errors(t)<-12
        errors(t)=errors(t)+24;
    end
end

% merrors = mean(errors);

allphases(11,1)= allphases(11,1)-24;
allphases(12,1)= allphases(12,1)-24;
allphases(13,1)= allphases(13,1)-24;
% 
allphases(11,3)= allphases(11,3)-24;
allphases(12,3)= allphases(12,3)-24;
allphases(11,9)= allphases(11,9)-24;
allphases(14,1)= allphases(14,1)-24;
%  

% what folows plots the array of figures that show for each gene and
% individual the residual (x-axis) against the Cosinor phase of the gene
% (y-axis)
allphases(3,1)= allphases(3,1)-24;
cm=hsv(10);
figure()
hold on
for k=1:num_probes; % running through probes
    subplot(4,4,k)
    plot(mean(errors),allphases(k,:),'*','Markersize',8,'Linewidth',1.5)
    lsline
    for indiv=1:10;
        hold on
        plot(mean(errors(:,indiv)),allphases(k,indiv),'*','Color',cm(indiv,:),'Markersize',8,'Linewidth',1.5)
    end

    [f gof]=fit(mean(errors)',allphases(k,1:10)','poly1');

    str=sprintf('f(x)= %.1f x+ %.1f R2= %.2f', f.p1 ,f.p2,gof.rsquare);
    hold on
    title({probenames{k}; str})
    plot([0 0],[min(allphases(k,:))-1 max(allphases(k,:))+1],'k-')
    plot([-4 3],[f.p2 f.p2],'k-')
    axis([-4 3 min(allphases(k,:))-1 max(allphases(k,:))+1])

    xlabel('Residual')
    ylabel('Cosinor phase')
    box on
    grid on
end
% end of plot