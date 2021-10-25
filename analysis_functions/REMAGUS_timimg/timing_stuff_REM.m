% this is the script for creating the plot that compares the REMAGUS data
% which has timing with the estimated timing coming from TT. 

Th_e = human_TT(esetfrmaBC,[15],-10);
Th=Th_e;
%%
load('times from Zeitg analysis')
%% tms has real time in first col, TT times in second and Theta in third, fourth is index of patient id and fifth is the time given by ZeitZeiger
[ r,tms ] = insert_real_times( remaguswithclock, Th.Tclock, Th.D_Thetas, 24*Ztimes_REM );
%% tms2 is the tms for all with a second peak
tms2 = tms(find(tms(:,2)~=-1),:);
%% this put the timing of the second peak into column 6 and the timing ofthe first peak if there is no second peak
ths=[];
for i=1:size(tms2,1)
    tms2(i,6)=mod(8+24*Th.second_Ts(floor(tms2(i,4))),24);
    if Th.second_Ts(floor(tms2(i,4))) == -1
        tms2(i,6)=tms2(i,2);
        ths = [ths Th.D_Thetas(floor(tms2(i,4)))];
    end
end
%% plot TT timing agaist real time
figure
scatter(tms(:,1),tms(:,2))
%% if TT timimg is too early replace it by timing of second peak in stms
stms=tms2;
for i=1:size(tms2,1)
    if tms2(i,2)<10
        stms(i,2)=tms2(i,6);
    end
end
%% 
figure
plot([0 24],[0 24],'k')
hold on
scatter(tms2(:,1),tms2(:,2),'b')
scatter(stms(:,1),stms(:,2),'r')
xlim([0 24])
ylim([0 24])
%%

%%
G=find(stms(:,2)<10);
inds = floor(stms(G,4));
figure
demoC(Th.Likelis(:,G))
%% for comparison with Zeitge Ztimes
figure
plot([0 24],[0 24],'k')
hold on
scatter(tms(:,2),tms(:,5),60,'b')
scatter(stms(:,2),stms(:,5),60,'r')
xlim([0 24])
ylim([0 24])
%%
load('timesTrain_1_15:01:17.mat')
%% This is for comparing Zeit and TT timing
% Run the section below first
% 0 is 8pm
cm=hsv(10);
rpred=reshape(out.timePred,6,10);
ract=reshape(out.timeObs,6,10);
hold on
for in=1:10
       scatter(mod(24*ract(:,in)-8,24)+1,mod(24*rpred(:,in)-8,24),60,cm(in,:),'filled','d');
end
grid on
box on
xlim([0 24])
ylim([0 24])
%plot([0 24],[0 24],'k')
%%
figure
cm=hsv(10);
plot([0 24],[0 24],'k')
hold on
tims=(Th_out.D_Ts*24/385);
for in=1:10
       scatter(0:4:20,tims(:,in),60,cm(in,:),'filled');
end
%%
[rho,pval] = corr(sThetas',lml(I)','Type','Pearson');
%%
er=abs(stms(:,1)-stms(:,2));
figure
histogram(er)
mer=sum(er)/length(er)
J=find(stms(:,2)>10);
merJ=sum(er(J))/length(er(J))
%% load likelis for Zeitzeiger
cd '/Users/davidrand/Dropbox/_Documents/_Work/Computing/R/TT stuff/v1/zeitzeiger/z1/ZanalysisREM'
load('likeliTest_1_11:58:42.mat')
likelis=out';
cd '/Users/davidrand/Dropbox/_Documents/_Work/Computing/Matlab_Stuff/time-teller_repository/DAR TT v2/TT_20-12/'
%% calculates Thetas for ZeitZeiger
Th_Z = Theta_calculator(likelis)
%% plots figure where Z Thetas and TT Thetas are plotted for each sample when the TT thetas are ordered,
x=(1:226);
[sThetas, I]=sort(Th.D_Thetas);
% get Z Thetas
ZThetas=Th_Z.Thetas;
plot(sThetas,'ro')
hold on
plot(ZThetas(I),'bo')
xlim([0 226])
ylim([0 1])
%%
[rho,pval] = corr(ZThetas',Th.D_Thetas','Type','Pearson');
%%
tbl = table(ZThetas',(1:226)')
lm = fitlm(tbl,'linear')
%%
lml=Th_Z.Log_max_likli;
figure
plot(sThetas,'ro')
hold on
plot(lml(I),'bo')
xlim([0 226])
%ylim([0 1])
%%
[rho,pval] = corr(sThetas',lml(I)','Type','Pearson');
%%
cd '/Users/davidrand/Dropbox/_Documents/_Work/Computing/R/TT stuff/v1/zeitzeiger/z1/ZanalysisREM'
load('times from Zeitg analysis REM')
%likelis=out';
cd '/Users/davidrand/Dropbox/_Documents/_Work/Computing/Matlab_Stuff/time-teller_repository/DAR TT v2/TT_20-12/'
%% plot of Z times against TT times for REM
% this figure is currently S14
% unmodified times
% is second_Ts is NaN curve is flat, if -1 there is only one peak
times_TT=mod(8+24*Th_e.Ts,24);
times_Z=24*Ztimes_REM;
times = zeros(length(times_TT),3);
times(:,1)=times_TT;
times(:,2)=times_Z;
times(:,3)=mod(8+24*Th_e.second_Ts,24);
% remove data with no second peak
KK=find(ismissing(times(:,3))~=0);
times2=times;
times2(KK,:)=[];
KK8=find(times2(:,3)==8);
times3=times2;
times3(KK8,:)=[];
figure
hold on
scatter(times3(:,1),times3(:,2),'b')
plot([0 24],[0 24],'k')
xlim([0 24])
ylim([0 24])
box on
grid on
KK=find(times3(:,1)<10);
times4=times3;
times4(KK,1)=times4(KK,3);
figure
scatter(times4(:,1),times4(:,2),'b')
hold on
scatter(times4(KK,1),times4(KK,2),'r')
% add in unimodal peaks
KK=find(times(:,3)==8)
scatter(times(KK,1),times(KK,2),'b')
plot([0 24],[0 24],'k')
xlim([0 24])
ylim([0 24])
box on
grid on
%%
xx=times3(KK,1);yy=times3(KK,2);
[rho,pval] = corr(xx,yy,'Type','Spearman')
%%
K=find(Th_e.Ts~=-1);
times_TT=times_TT(K);
times_Z = times_Z(K);
figure
hold on
scatter(times_TT,times_Z,'b')
%% move using second peaks
hold on
second_timesTT=Th_e.second_Ts(K);
J=find(and(times_TT<10,second_timesTT~=-1));
times_TT_altered = times_TT;
times_TT_altered(J) = mod(8+24*Th_e.second_Ts(J),24);
scatter(times_TT(J),times_Z(J),'r')
scatter(times_TT_altered(J),times_Z(J),'g')
xlim([0 24])
ylim([0 24])
figure
K=find(and(times_TT_altered~=8,times_TT_altered~=NaN));
scatter(times_TT_altered(K),times_Z(K),'k')
xlim([0 24])
ylim([0 24])
xx=times_TT_altered(K)';
KK=find(ismissing(xx)==0);
yy = times_Z(K);
xx = xx(KK);
yy = yy(KK);
JJ=find(yy>9);
xx=xx(JJ);
yy=yy(JJ);
%figure
scatter(xx,yy,'r')
[rho,pval] = corr(xx,yy,'Type','Pearson')
xlim([0 24])
ylim([0 24])
%% maen abs error of Ztimes
mean(abs(tms(:,1)-tms(:,5)))
%% mean abs error of TTtimes
mean(abs(tms(:,1)-tms(:,2)))
%% mean abs error of TTtimes after alteration
mean(abs(stms(:,1)-stms(:,2)))
%%
