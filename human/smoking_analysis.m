% Smoking analysis
% difference in smoking and non-smoking Thetas is significant
%%
% load dataset
load('esetfrmaSMOM_data')
%Smoker data: 1:40 non-smoker, 41:79 smoker
%%
% calculate TT stuff
Th_smoke = human_TT(esetfrmaSMOM_data,[11 15],-12);
%%
% make boxplots
group4=[zeros(1,40) ones(1,39)];
boxplot(Th_smoke.D_Thetas,group4);
%%
% ttest
% cancer vs normal
[h,p,ci,stats] = ttest2(Th_smoke.D_Thetas(1:40),Th_smoke.D_Thetas(41:79),'Vartype','unequal')% cancer & Dysplasic against normal
% save the Th structure
%save('Feng Th Structure','Th_Feng')
%%
% gives an interesting plot
demoC(Th_smoke.Likelis)
%%