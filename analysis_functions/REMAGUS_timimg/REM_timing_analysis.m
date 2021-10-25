function s = REM_timing_analysis(Ts,Thetas,Ztimes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load('remagus_factor_inc_clockfunc.mat');

rwc=remaguswithclock;
[rwc2,rts] = insert_real_times( rwc,Ts,Thetas,Ztimes );
% This produces a table with the real times inserted and also outputs an
% array with the real times in the first column, the estimate times in
% the second and the Thetas in the third.
K=find(rts(:,3)<0.05);
rts2=rts(K,:);
figure
scatter(rts(:,1),rts(:,2),'r');
hold on
scatter(rts2(:,1),rts2(:,2),'b');
plot([0 24],[0 24],'k')
hold off
xlim([0 24]);
ylim([0 24]);

K2=find(rts2(:,2)<10);
K=find(rts(:,2)<10);
[p,h,stats] = ranksum(rts(K,2),rts2(K2,2),'tail','left');
s.p=p;
s.allThetas_goodtime=find(rts(:,2)>=10);
s.allThetas_badtime=K;
s.goodThetas_goodtime=find(rts2(:,2)>=10);
s.goodThetas_badtime=K2;
end

