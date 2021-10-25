%%
Th_out = leave_1_out_Bjarn([],-10);
%%
load('Zhang Th Structure 5')
Th=Th_Zhang_5;
%%
eight_tissues=[1 2 4 6 8 9 10 11];
load('zhang11_12organs_data')
%%
organnames={'Adrenal','Aorta','BrownFat','Heart','Kidney','Liver','Lung','SkelMus'};
times = {'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty','twentyone','twentytwo','twentythree','twentyfour'};

for tms = 1:24; % times
    counter=0;
    for j=[1 2 4 6 8 9 10 11]; % 8 tissues
        counter = counter+1; % allows for possibility that some js are missing
        clockstr.(times{tms})(:,counter) = zhang11_12organs(:,24*(j-1)+tms);
        % clock11hogenstr.one is a 11x8 array, 11 probe expressions times 8 timepoints
    end
end
%%
data = zhang11_12organs;
data(g,a)=data(g,[(a-1)*6+1:6*a, (a-1)*6+1])
%%

datastr = make_clock_OM_str(zhang11_12organs,rhythmic_probeset,all_probe_names,instances_list,timeset,data_type)

%%
% use cosinor to get the phases of the probes for each individual
clear q a p allphases maxf minf antiphase
%%
% q probe number
% a instances(individuals, tissues)
num_probes=11;
odds=1:2:23;
evens=2:2:24;
for q=1:num_probes; % runs over probesnum_probes
    for a = 1:8; % runs over individuals
        data=[];
        for i = 1:24
            data = [data, clockstr.(times{i})(q,a)];
        end
        t4data=data;
        t4data(evens)=data(13:24);
        t4data(odds)=data(1:12);
        [p(q,a) allphases(q,a),maxf(q,a),antiphase(q,a), minf(q,a)]= cosinor(linspace(0,1,24),t4data,2*pi,0.05);
        %[p(q,a) allphases(q,a),maxf(q,a),antiphase(q,a), minf(q,a)]= cosinor(linspace(0,1,7),clock_OM(q,[(a-1)*6+1:6*a, (a-1)*6+1]),2*pi,0.05);

        if allphases(q,a)>24
            allphases(q,a)=allphases(q,a)-24;
        end
    end
end
%%
ac=repmat(actuals,1,8);
errors=ac-24*Th.Ts;
errmod24=mod(errors,24);
errors=min(errmod24,abs(errmod24-24));
%%
mns=mean(allphases')';
mns2=repmat(mns,1,8);
allphases-mns2