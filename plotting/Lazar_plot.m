%% Lazar_plot
% this analyses and plots things assocoiated with the Fang/Lazar data
load('esetfrmaLazar')
Th = mouse_TT( esetfrmaLazar,[6],24,-5);
group=[ones(1,5) 2*ones(1,5)];
[p,h,stats] = ranksum(Th.D_Thetas(1:5),Th.D_Thetas(6:10),'Tail','left')
%%
figure;boxplot(Th.D_Thetas(:),group)
%%
x=24*(1:193)/193;
figure
plot(x+18,log(Th.Likelis(:,1:5)),'b')
hold on
plot(x+18,log(Th.Likelis(:,6:10)),'r')
xlim([15 45])
%%
figure
plot(x+18,Th.Likelis(:,1:5),'b')