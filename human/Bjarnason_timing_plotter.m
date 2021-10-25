function Bjarnason_timing_plotter(Th)
% plots the timing scatter plot of actual vs estimated
num_finetimes = size(Th.Likelis_all,2);


actuals = linspace(0,20,6)';
actuals = [8 12 16 20 0 4]';
figure()
fa = 24/num_finetimes;
Ts=Th.D_Ts;
Ts=mod(fa*Ts+8,24);
hold on
clrs = distinguishable_colors(10,'w');
for individual=1:10;
    plot(actuals,Ts(:,individual),'o','Color',clrs(individual,:));
end
% plot(actuals, Ts(:,1),'bo','Markersize',8)
% plot(actuals, Ts(:,2),'b*','Markersize',8)
% plot(actuals, Ts(:,6:10),'m*','Markersize',8)
% plot(actuals, Ts(:,3:5),'b*','Markersize',8)
plot([-1 25], [-1 25],'k-')
plot([-1 1], [23 25],'k-')
plot([23 25], [-1 1],'k-')
grid on
box on
axis([ -1 25 -1 25])
xlabel('Real Time')
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
ylabel('Estimated Time')
set(gca,'YTick' ,0:4:24)
set(gca,'YTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
legend( 'Phase shifted Male','Male','Female')

    % this will plot the Bjarnason data into the top twp principal
    % components for each timepoint

% plots the likelihood functions
Mall=Th.Likelis_all;
%a=a2/num_finetimes*24;
cm2=clrs;%jet(10);
figure()
hold on
for individual=1:10;
    for tms =1:6;
        hold on
        subplot(2,3,tms)
        M2 = reshape(Mall(individual,:,tms),num_finetimes,1);
        p=area(linspace(0,24,num_finetimes),M2,'Facecolor',cm2(individual,:));
        alpha(p,0.2)
    end
end
legend('M15','M11','M12','M11','M9','F18','F13','F14','F5','F6')
for tms =1:6;
	hold on
    subplot(2,3,tms)
    plot([(tms-1)*4,(tms-1)*4],[0,6],'r-','Linewidth',3)
    axis([0 24 0 max(max(Mall(:,:,tms)))])
    set(gca,'XTick' ,0:4:24)
    xtickangle(90)
    set(gca,'XTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
    xlabel('Time')
    ylabel('Likelihood')
end
subplot(2,3,1)
plot([(7-1)*4,(7-1)*4],[0,6],'r-','Linewidth',3)

% plots a histogram of Thetas
max(max(Th.D_Thetas))
figure()
hold on
m=max(Th.D_Thetas(:))
histogram(Th.D_Thetas(:),0:0.005:m)
box on
grid on
ylabel('Count')
xlabel('\Theta')

% plos more likelihood functions: set up verticals
%a=a2/num_finetimes*24;
cm2=hsv(6);
figure()
hold on
individual=1;
for tms =1:6;
    hold on
    p=area(linspace(0,24,num_finetimes),reshape(Mall(individual,:,tms),num_finetimes,1),'Facecolor',cm2(tms,:));
    alpha(p,0.2)
    plot([(tms-1)*4,(tms-1)*4],[0,6],'-','Linewidth',3,'Color',cm2(tms,:))
end
plot([(7-1)*4,(7-1)*4],[0,6],'-','Linewidth',3,'Color',cm2(1,:))
xtickangle(90)
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
xlabel('Time')
ylabel('Likelihood')
%axis([0 24 0 5.7])
box on

% now plot likelihoods
for tms =1:6;
    hold on
    subplot(2,3,tms)
    plot([(tms-1)*4,(tms-1)*4],[0,6],'r-','Linewidth',3)
    axis([0 24 0 max(max(Mall(:,:,tms)))])
    set(gca,'XTick' ,0:4:24)
    xtickangle(90)
    set(gca,'XTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
    xlabel('Time')
    ylabel('Likelihood')
end
subplot(2,3,1)
plot([(7-1)*4,(7-1)*4],[0,6],'r-','Linewidth',3)
return
