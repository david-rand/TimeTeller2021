%% Oster timing stuff
Th = mouse_TT( esetfrmaOstercore,[],24,-5 );
%%
Th=Th12;
%% This plots the timing scatter figure for the Oster data
figure
a=[1 1 1 7 7 7 13 13 13 19 19 19];
actuals = [a a a a];
x=mod(actuals+12,24);y=mod(8+24*Th.Ts,24);
plot(x(1:12), y(1:12),'bo','MarkerSize',10,'LineWidth',1)
hold on
plot([0 24],[0 24],'k')
plot(x(25:36), y(25:36),'or','MarkerSize',10,'LineWidth',1)
inds=[34 35 36];
yy=mod(8+24*Th.second_Ts(inds),24);
xx=x(inds)
plot(xx,yy,'om','MarkerSize',10,'LineWidth',1);
xlim([0 24])
ylim([0 24])
box on
grid on
%% Thsi does the boxplot for the Oster data
figure
prob=Th.D_Thetas;
group3=[ones(1,12) 2*ones(1,12) 3*ones(1,12) 4*ones(1,12) ];
plotSpread([prob(1:12);prob(13:24);prob(25:36);prob(37:48)]','distributionIdx',[ones(1,12), 2*ones(1,12),3*ones(1,12), 4*ones(1,12)],'distributionMarkers',{'*','*','*','*'},'distributionColors',{'k','k','k','k'})
boxplot(prob,group3)
%% Ranksum analysis
[p,h,stats] = ranksum(Th.D_Thetas(1:12),Th.D_Thetas(13:24),'Tail','left')
%% This calculates the absolute errors
e=abs(y-x);e1=e(1:12);e2=e(25:36);
disp('abs error')
disp((sum(e1)+sum(e2))/length(x));
yyy=y;
yyy(inds)=yy;
e=abs(yyy-x);e1=e(1:12);e2=e(25:36);
disp('abs error after adjustment')
disp((sum(e1)+sum(e2))/length(x));
%%