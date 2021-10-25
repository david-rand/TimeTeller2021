%% loads
load('esetfrmaBC Th Structure 12.mat')
load('Feng Th Structure excl15.mat')
load('smoking Th Structure 12 excl11&15.mat')
load('frmaOM Th Structure 12.mat')
%%
% put Thetas together
Ths.B=Th_frmaOM.D_Thetas; 
Ths.S=Th_smoke.D_Thetas(1:40); 
Ths.FN=Th_Feng.D_Thetas(168:212);
Ths.OC=Th_Feng.D_Thetas(1:167); 
Ths.DYS=Th_Feng.D_Thetas(213:229);
Ths.REM=Th_esetfrmaBC.D_Thetas;
% Thetas = [Th_frmaOM.D_Thetas Th_smoke.D_Thetas(1:40) Th_Feng.D_Thetas(168:212) Th_Feng.D_Thetas(1:167) Th_Feng.D_Thetas(213:229) Th_esetfrmaBC.D_Thetas];
Thetas = [Ths.B Ths.S Ths.FN Ths.OC Ths.REM];
%%
% make groups
group_Bjarn=6*ones(1,length(Th_esetfrmaBC.D_Thetas));
group_feng_dys=5*ones(1,230-213);
group_feng_oc=4*ones(1,167);
group_feng = 3*ones(1,213-168);
group_smoke=2*ones(1,40);
group_frma=ones(1,60);
%group = [group_frma group_smoke group_feng group_feng_oc group_feng_dys group_Bjarn];
group = [group_frma group_smoke group_feng group_feng_oc group_Bjarn];
%%
% Plot
boxplot(Thetas,group)
%%
nma_small={'B','S','FN'};
nma_big={'OC','DYS','REM'};

% test means
for i=1:3 % bigger
    for j=1:3
        [h,p] = ttest2(Ths.(nma_big{i}),Ths.(nma_small{j}),'Vartype','unequal','Tail','right'); % tests x greater than y
        pv.(nma_big{i}).(nma_small{j})=p;
    end
end
%%

