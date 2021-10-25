% Feng analysis
%%
% load dataset
load('Fengfrma')
%items 1:167 are o cancer, 213:229 dysplasic and 168:212 normal
%%
% calculate TT stuff
Th_Feng = human_TT( Fengfrma,[],-12)
%%
% make boxplots
figure
group = [3*ones(1,167) ones(1,45),2*ones(1,17)];
boxplot(Th.D_Thetas,group);
%%
% ttest
% cancer vs normal
[h1,p1,ci1,stats1] = ttest2(Th_Feng.D_Thetas(1:167),Th_Feng.D_Thetas(168:212),'Vartype','unequal','Tail','right');
[p2,h2,stats2] = ranksum(Th_Feng.D_Thetas(1:167),Th_Feng.D_Thetas(168:212),'Tail','right');
% cancer & Dysplasic against normal
[h3,p3,ci3,stats3] = ttest2([Th_Feng.D_Thetas(1:167) Th_Feng.D_Thetas(213:229)],Th_Feng.D_Thetas(168:212),'Vartype','unequal','Tail','right');
[p4,h4,stats4] = ranksum([Th_Feng.D_Thetas(1:167) Th_Feng.D_Thetas(213:229)],Th_Feng.D_Thetas(168:212),'Tail','right');
% normal vs dysplasic
[h5,p5,ci5,stats5] = ttest2(Th_Feng.D_Thetas(213:229),Th_Feng.D_Thetas(168:212),'Vartype','unequal','Tail','right');
[p6,h6,stats6] = ranksum([Th_Feng.D_Thetas(213:229)],Th_Feng.D_Thetas(168:212),'Tail','right');
%dys vs cancer
[p7,h7,stats7] = ranksum([Th_Feng.D_Thetas(1:167) ],Th_Feng.D_Thetas(213:229),'Tail','left');

%% save the Th structure
save('Feng Th Structure','Th_Feng')
%
%%
% gives an interesting plot
demoC(Th_Feng.Likelis)
%%
figure
hold on
tmsfeng=mod(8+24*Th.Ts,24);
scatter(tmsfeng(1:167),Th.D_Thetas(1:167),60,'r')
scatter(tmsfeng(168:212),Th.D_Thetas(168:212),60,'b')
scatter(tmsfeng(213:229),Th.D_Thetas(213:229),60,'k')
JJ = find(tmsfeng<7);
KK = find(tmsfeng>=7);
tms2feng = mod(8+24*Th.second_Ts,24);
scatter(tms2feng(JJ),Th.D_Thetas(JJ),60,'g')
%xlim([0 24])
%%
Y=[Th_Feng.D_Thetas(168:212)',Th_Feng.D_Thetas(213:229)',Th_Feng.D_Thetas(1:167)'];
figure
[h,L,MX,MED]=violin(Y);

