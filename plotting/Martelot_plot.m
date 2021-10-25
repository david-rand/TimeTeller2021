function Martelot_plot(Th)

%timing scatter plot

%% this is the order the data is in


numfinetimes = size(Th.Likelis,1);

actuals = [ 2 6 10 14 18 26 22];
% 0 is 8pm so to get the actual times above we use
rt = mod(-6+Th.D_Ts*24/numfinetimes,24);

scatter(actuals,rt,'*','r')
hold on
plot([-1 49], [-1 49],'k-')
plot([-1 25], [23 49],'k-')
plot([23 49], [-1 25],'k-')
plot([-1 1], [47 49],'k-')
plot( [47 49],[-1 1],'k-')
axis([0 27 0 27])
xlabel('Real Time')
ylabel('Estimated Time')
box on
grid on


% plot likelihood and log likelihood curves
xx=24*(1:193)/193;
rts = mod(-6+xx,24);
[srts,I]=sort(rts);
ll=Th.Likelis./max(Th.Likelis);
figure
plot(srts,Th.Likelis(I,:))
xlim([0 24])
figure
plot(srts,ll(I,:))
xlim([0 24])

% plot boxplot for Thetas
figure()
boxplot(Th.D_Thetas)
ylabel('Theta')
set(gca,'XTickLabel',{'WT'})

% This plots the Zhang et al timeseries and compares the LeMartelot data
% figure()
% for i = 1:11;
%     subplot(4,3,i)
%     hold on
%           for j=[1 2 4 6 8 9 10 11 12];       
%  plot(18:2:64,hogen11(i,24*(j-1)+1:24*j),'k-');
%           end
% plot(2+24,clock11Martletot(i,1),'bo') 
% plot(6+24,clock11Martletot(i,2),'ro') 
% plot(10+24,clock11Martletot(i,3),'mo')
% plot(14+24,clock11Martletot(i,4),'go') 
% plot(18+24,clock11Martletot(i,5),'b*') 
% plot(22+24,clock11Martletot(i,6),'r*') 
% plot(26+24,clock11Martletot(i,7),'m*')
% axis([18 66 5 12])
% title(Inx11names(i))
% xlabel('Time')
% ylabel('Expression')
end