function make_REMAGUS_csv(Thetas, Ts, filename)

% This writes a csv file with the REMAGUS data table in it.
% Inputs: Thetas: the Theta values for the given dataset, usually from
% human_TT
% Ts: the corresponding times
% filename: the name of the output file
% USAGE: e.g. make_REMAGUS_csv(D_thetas,Ts,'rem_1.csv')

load('remagus_factor_inc_clockfunc.mat')
remaguswithclock.ClockFunction2=Thetas';
% lists Thetas for the 226 samples
remaguswithclock.estimatedtime=Ts';
% lists the times for the 226 samples
% 
% m=max(Thetas);
% 
% figure(14)
% hold on
% histogram(prob,0:0.05:m)

% all the stuff here is concerned with assigning a numeric value to
% different pronostic factors
% record tumor size into structure
for i=1:226;  
    if remaguswithclock.TumourSize(i)=='2';
        remaguswithclock.TumourSize3(i)=0;
    else
        if remaguswithclock.TumourSize(i)=='3';
            remaguswithclock.TumourSize3(i)=1;
        else
            if remaguswithclock.TumourSize(i)=='4';
                remaguswithclock.TumourSize3(i)=2;
            else
                remaguswithclock.TumourSize3(i)=0;
            end
        end
    end
end

for i=1:226; % this is better done using switch cases
    if remaguswithclock.NodalStatus2(i)=='0';
        remaguswithclock.NodalStatus3(i)=0;
    else
        if remaguswithclock.NodalStatus2(i)=='1';
            remaguswithclock.NodalStatus3(i)=1;
        else
            if remaguswithclock.NodalStatus2(i)=='2';
                remaguswithclock.NodalStatus3(i)=2;
            else
                if remaguswithclock.NodalStatus2(i)=='3';
                    remaguswithclock.NodalStatus3(i)=3;
                else
                    remaguswithclock.NodalStatus3(i)=4;
                end
            end
        end
    end
    
    if remaguswithclock.HER2status(i)=='HER2 pos';
        remaguswithclock.HER2status2(i)=1;
    else
        if  remaguswithclock.HER2status(i)=='HER2 neg';
            remaguswithclock.HER2status2(i)=0;
        else
       remaguswithclock.HER2status3(i)='Nan';
        end
    end

    if remaguswithclock.ERstatus(i)=='ER pos';
        remaguswithclock.ERstatus2(i)=1;
    else
        if  remaguswithclock.ERstatus(i)=='ER neg';
            remaguswithclock.ERstatus2(i)=0;
        else
       remaguswithclock.ERstatus2(i)='Nan';
        end
    end

    if remaguswithclock.PRstatus(i)=='PR pos';
        remaguswithclock.PRstatus2(i)=1;
    end
    if  remaguswithclock.PRstatus(i)=='PR neg';
        remaguswithclock.PRstatus2(i)=0;
    end
    if remaguswithclock.PRstatus(i)=='#N/A';
        remaguswithclock.PRstatus2(i)=2;
    end

    if remaguswithclock.pCR2(i)=='pCR';
        remaguswithclock.pCR3(i)=1;
    end
    if  remaguswithclock.pCR(i)=='No pCR';
        remaguswithclock.pCR2(i)=0;
    end
    if  remaguswithclock.pCR(i)=='#N/A';
        remaguswithclock.pCR2(i)=2;
    end
end

for i=1:226;
    if  remaguswithclock.HER2status2(i)==0 && remaguswithclock.ERstatus2(i)==0 && remaguswithclock.PRstatus2(i)==0;
        remaguswithclock.TriNeg2(i)=1;
    elseif   remaguswithclock.HER2status2(i)==2 | remaguswithclock.ERstatus2(i)==2 | remaguswithclock.PRstatus2(i)==2;
        remaguswithclock.TriNeg2(i)=2;
    else
         remaguswithclock.TriNeg2(i)=0;
    end
end

for i=1:226;
	if  remaguswithclock.HER2status2(i)==1 && remaguswithclock.ERstatus2(i)==1 && remaguswithclock.PRstatus2(i)==0;
        remaguswithclock.PRNeg2(i)=1;
    elseif   remaguswithclock.HER2status2(i)==2 | remaguswithclock.ERstatus2(i)==2 | remaguswithclock.PRstatus2(i)==2;
        remaguswithclock.PRNeg2(i)=2;
    else
        remaguswithclock.PRNeg2(i)=0;
    end
   %%%%%%%
    if  remaguswithclock.HER2status2(i)==1 && remaguswithclock.ERstatus2(i)==0 && remaguswithclock.PRstatus2(i)==1;
        remaguswithclock.ERNeg2(i)=1;
    elseif   remaguswithclock.HER2status2(i)==2 | remaguswithclock.ERstatus2(i)==2 | remaguswithclock.PRstatus2(i)==2;
        remaguswithclock.ERNeg2(i)=2;
    else
     	remaguswithclock.ERNeg2(i)=0;
  	end
    %%%%%%%
	if  remaguswithclock.HER2status2(i)==0 && remaguswithclock.ERstatus2(i)==1 && remaguswithclock.PRstatus2(i)==1;
        remaguswithclock.HER2Neg2(i)=1;
    elseif   remaguswithclock.HER2status2(i)==2 | remaguswithclock.ERstatus2(i)==2 | remaguswithclock.PRstatus2(i)==2;
        remaguswithclock.HER2Neg2(i)=2;
    else
        remaguswithclock.HER2Neg2(i)=0;
    end
end

for i=1:226;
    if  remaguswithclock.HER2status2(i)==0 && remaguswithclock.ERstatus2(i)==0 && remaguswithclock.PRstatus2(i)==1;
        remaguswithclock.PRPos2(i)=1;
    elseif   remaguswithclock.HER2status2(i)==2 | remaguswithclock.ERstatus2(i)==2 | remaguswithclock.PRstatus2(i)==2;
    	remaguswithclock.PRPos2(i)=2;
    else
       	remaguswithclock.PRPos2(i)=0;
    end
   %%%%%%%
   	if  remaguswithclock.HER2status2(i)==0 && remaguswithclock.ERstatus2(i)==1 && remaguswithclock.PRstatus2(i)==0;
        remaguswithclock.ERPos2(i)=1;
    elseif   remaguswithclock.HER2status2(i)==2 | remaguswithclock.ERstatus2(i)==2 | remaguswithclock.PRstatus2(i)==2;
        remaguswithclock.ERPos2(i)=2;
    else
       	remaguswithclock.ERPos2(i)=0;
   	end
    %%%%%%%
   	if  remaguswithclock.HER2status2(i)==1 && remaguswithclock.ERstatus2(i)==0 && remaguswithclock.PRstatus2(i)==1;
        remaguswithclock.HER2Pos2(i)=1;
    elseif   remaguswithclock.HER2status2(i)==2 | remaguswithclock.ERstatus2(i)==2 | remaguswithclock.PRstatus2(i)==2;
   	remaguswithclock.HER2Neg2(i)=2;
    else
       	remaguswithclock.HER2Pos2(i)=0;
    end
end



for i=1:226;
    if  remaguswithclock.pCR2(i)==0 && remaguswithclock.TriNeg2(i)==0 ;
        remaguswithclock.PCRTN(i)=1;
    elseif   remaguswithclock.pCR2(i)==0 && remaguswithclock.TriNeg2(i)==1;
   	remaguswithclock.PCRTN(i)=2;
  	elseif   remaguswithclock.pCR2(i)==1 && remaguswithclock.TriNeg2(i)==0;
  	remaguswithclock.PCRTN(i)=3;
    else
      	remaguswithclock.PCRTN(i)=4;
    end
end
writetable(remaguswithclock,filename);
%return

prob=Thetas;

for figloop = 1
    % ER status box
    figure()
    hold on
    count_estro0 = sum(remaguswithclock.ERstatus2==0);
    count_estro1 = sum(remaguswithclock.ERstatus2==1);
    p_ER =ranksum(prob(remaguswithclock.ERstatus2==0),prob(remaguswithclock.ERstatus2==1))
    subplot(4,2,1)
    plotSpread([prob'],'distributionIdx',[remaguswithclock.ERstatus2'+1],'distributionMarkers',{'.','.'},'distributionColors',{'k','k'})
    boxplot([prob],[remaguswithclock.ERstatus2])
    title({ sprintf('Estrogen Receptor p=%0.3f' , p_ER)})
    set(gca,'XTickLabel',{sprintf('Negative n=%d' , count_estro0),sprintf('Positive n=%d' , count_estro1)})
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    axis( [0.5 2.5 -inf inf])
    box on
    ylabel('\Theta')

    % PR status box
     p_PR =ranksum(prob(remaguswithclock.PRstatus2==0),prob(remaguswithclock.PRstatus2==1))
    count_progestro0 = sum(remaguswithclock.PRstatus2==0);
    count_progestro1 = sum(remaguswithclock.PRstatus2==1);
    subplot(4,2,2)
    plotSpread([prob'],'distributionIdx',[remaguswithclock.PRstatus2'+1],'distributionMarkers',{'.','.','.'},'distributionColors',{'k','k','k'})
    boxplot(prob,[remaguswithclock.PRstatus2])
    set(gca,'XTickLabel',{sprintf('Negative n=%d' , count_progestro0),sprintf('Positive n=%d' , count_progestro1)})
    axis( [0.5 2.5 -inf inf])
    title({ sprintf('Progesterone Receptor p=%0.3f' , p_PR)})
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    ylabel('\Theta')

    % Triple neg box status box
    count_trineg0 = sum(remaguswithclock.TriNeg2==0);
    count_trineg1 = sum(remaguswithclock.TriNeg2==1);
    p_TN =ranksum(prob(remaguswithclock.TriNeg2==0),prob(remaguswithclock.TriNeg2==1))
    subplot(4,2,4)
    plotSpread([prob'],'distributionIdx',[remaguswithclock.TriNeg2'+1],'distributionMarkers',{'.','.','.'},'distributionColors',{'k','k','k'})
    boxplot([prob],[remaguswithclock.TriNeg2])
    title({ sprintf('Triple Negative p=%0.4f' , p_TN)})
    set(gca,'XTickLabel',{sprintf('No n=%d' , count_trineg0),sprintf('Yes n=%d' , count_trineg1)})
    axis( [0.5 2.5 -inf inf])
    box on
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    ylabel('\Theta')

    % HER2 status box
    count_Her2pos = sum(remaguswithclock.HER2status=='HER2 pos');
    count_Her2neg = sum(remaguswithclock.HER2status=='HER2 neg');
    p_Her2 =ranksum(prob(remaguswithclock.HER2status=='HER2 pos'),prob(remaguswithclock.HER2status=='HER2 neg'))
    subplot(4,2,3)
    plotSpread([prob'],'distributionIdx',[remaguswithclock.HER2status2'+1],'distributionMarkers',{'.','.'},'distributionColors',{'k','k'})
    boxplot([prob],[remaguswithclock.HER2status2])
    title({ sprintf('HER2 p=%0.4f' , p_Her2)})
    set(gca,'XTickLabel',{sprintf('Negative n=%d' , count_Her2neg ),sprintf('Positive n=%d' , count_Her2pos)})
    axis( [0.5 2.5 -inf inf])
    box on
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    ylabel('\Theta')


    % now three Wilcoxon rank sum tests
    [p_gr1] =ranksum(prob(remaguswithclock.Grade==1),prob(remaguswithclock.Grade==2))
    [p_gr2] =ranksum(prob(remaguswithclock.Grade==1),prob(remaguswithclock.Grade==3))
    [p_gr3] =ranksum(prob(remaguswithclock.Grade==2),prob(remaguswithclock.Grade==3))

    % Grade box
    subplot(4,2,5)
    plotSpread([prob'],'distributionIdx',[remaguswithclock.Grade'],'distributionMarkers',{'.','.','.'},'distributionColors',{'k','k','k'})
    boxplot([prob],[remaguswithclock.Grade])
    set(gca,'XTickLabel',{'Grade 1 n=15','Grade 2 n=83','Grade 3 n=121'})
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    title({ sprintf('Grade p_{12}=%0.3f, p_{13}=%0.3f, p_{23}=%0.3f' , p_gr1,p_gr2,p_gr3)})
    grid on
    axis( [0.5 3.5 -inf inf])
    ylabel('\Theta')

    % Nodal status box
    subplot(4,2,6)
    count_node0 = sum(remaguswithclock.NodalStatus3==0);
    count_node1 = sum(remaguswithclock.NodalStatus3==1);
    count_node2 = sum(remaguswithclock.NodalStatus3==2);
    count_node3 = sum(remaguswithclock.NodalStatus3==3);
    p_node =ranksum([prob(remaguswithclock.NodalStatus3==0),prob(remaguswithclock.NodalStatus3==1)],[prob(remaguswithclock.NodalStatus3==2),prob(remaguswithclock.NodalStatus3==3)])
    plotSpread([prob'],'distributionIdx',[remaguswithclock.NodalStatus3'+1],'distributionMarkers',{'.','.','.','.','.'},'distributionColors',{'k','k','k','k','k'})
    boxplot([prob],[remaguswithclock.NodalStatus3])
    title({ sprintf('Nodal Status 0&1 vs 2&3 p=%0.4f' , p_node)})
    set(gca,'XTickLabel',{sprintf('0 n=%d' , count_node0 ),sprintf('1 n=%d' ,  count_node1 ),sprintf('2 n=%d' ,  count_node2 ),sprintf('3  n=%d' ,  count_node3 )})
    axis( [0.5 4.5 -inf inf])
    box on
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    ylabel('\Theta')

    % tumour size status box
    subplot(4,2,7)
    count_TumourSize2 = sum(remaguswithclock.TumourSize=='2');
    count_TumourSize3 = sum(remaguswithclock.TumourSize=='3');
    count_TumourSize4 = sum(remaguswithclock.TumourSize=='4');
    p_TumourSize =ranksum(prob(remaguswithclock.TumourSize=='2'),[prob(remaguswithclock.TumourSize=='3'),prob(remaguswithclock.TumourSize=='4')])
    plotSpread([prob'],'distributionIdx',[remaguswithclock.TumourSize3'+1],'distributionMarkers',{'.','.','.'},'distributionColors',{'k','k','k'})
    boxplot([prob],[remaguswithclock.TumourSize3])
    title({ sprintf('Tumour Size 2 vs 3&4 p=%0.4f' , p_TumourSize)})
    set(gca,'XTickLabel',{sprintf('2 n=%d' , count_TumourSize2 ),sprintf('3 n=%d' ,  count_TumourSize3 ),sprintf('4 n=%d' ,  count_TumourSize4 )})
    axis( [0.5 4.5 -inf inf])
    box on
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    ylabel('\Theta')

    %pCR status box
    count_PCR0 = sum(remaguswithclock.pCR2==0);
    count_PCR1 = sum(remaguswithclock.pCR2==1);
    p_pCR =ranksum(prob(remaguswithclock.pCR2==0),prob(remaguswithclock.pCR2==1))
    subplot(4,2,8)
    plotSpread([prob'],'distributionIdx',[remaguswithclock.pCR2'+1],'distributionMarkers',{'.','.','.'},'distributionColors',{'k','k','k'})
    boxplot([prob],[remaguswithclock.pCR2])
    title({ sprintf('pCR p=%0.4f' , p_pCR)})
    set(gca,'XTickLabel',{sprintf('No n=%d' , count_PCR0 ),sprintf('Yes n=%d' , count_PCR1)})
    axis( [0.5 2.5 -inf inf])
    box on
    set(findobj(gca,'type','line'),'linew',2)
    grid on
    ylabel('\Theta')
end
% end of figloop

%%

% makes figure plotting Theta distributions for combinations of pCR status
% and TN status

count_PCR0TN = sum(remaguswithclock.pCR2==0 & remaguswithclock.TriNeg2==1)
count_PCR0TN0 = sum(remaguswithclock.pCR2==0 & remaguswithclock.TriNeg2==0)
count_PCR1TN = sum(remaguswithclock.pCR2==1 & remaguswithclock.TriNeg2==1)
count_PCR1TN0 = sum(remaguswithclock.pCR2==1 & remaguswithclock.TriNeg2==0)

p_pCRTN =ranksum(prob(remaguswithclock.PCRTN==3),prob(remaguswithclock.PCRTN==4))

figure()
plotSpread([prob'],'distributionIdx',[remaguswithclock.PCRTN'],'distributionMarkers',{'.','.','.','.'},'distributionColors',{'k','k','k','k'})
boxplot([prob],[remaguswithclock.PCRTN])
title({ sprintf('pCR p=%0.4f' , p_pCR)})
set(gca,'XTickLabel',{sprintf('No pCR, not TN n=%d' , count_PCR0TN0 ),sprintf('No pCR, TN  n=%d' , count_PCR0TN ),sprintf('pCR, not TN  n=%d' , count_PCR1TN0  ),sprintf('pCR, TN  n=%d' , count_PCR1TN )})
%axis( [0.5 2.5 -inf inf])
box on
set(findobj(gca,'type','line'),'linew',2)
grid on
ylabel('\Theta')

load('heureswithalpha')
for i = 1:108
if heureswithalpha.Alpha(i)<0.1;
    heureswithalpha.col(i) = 1;
else
    heureswithalpha.col(i) =3;
end
end
figure;
alpcol=[0 0 1; 0 1 0; 1 0 0]
hold on
plot(-1, -1,'o','color',alpcol(1,:),'Markerfacecolor',alpcol(1,:))
plot(-3, -3,'o','color',alpcol(3,:),'Markerfacecolor',alpcol(3,:))

 for i = 1:108; plot(heureswithalpha.RealTime(i),heureswithalpha.VarName5(i),'o','color',alpcol(heureswithalpha.col(i),:),'Markerfacecolor',alpcol(heureswithalpha.col(i),:),'MarkerSize',9); hold on; end
hold on 
plot([0 24],[0 24],'k-')
grid on ;box on
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'00:00','04:00','8:00','12:00','16:00','20:00','00:00'})
set(gca,'YTick' ,0:4:24)
set(gca,'YTickLabel',{'00:00','04:00','8:00','12:00','16:00','20:00','00:00'})
xlabel('Real Time of sample')
ylabel('Estimated Time of sample')
axis([ 0 24 0 24])
legend({'\Theta<0.1', '0.1<\Theta'}) %'0.1<\alpha<0.155'

figure()
for i =1:108; hold on
plot(abs(heureswithalpha.VarName5(i)-heureswithalpha.RealTime(i)),heureswithalpha.Alpha(i),'b*')
end
xlabel('Absolute error of estimation (Hours)')
ylabel('\Theta')
box on
grid on