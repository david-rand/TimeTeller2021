function demoC(Mall)

% This takes a sxd array of likelihoods and plots there likelihood ratio curves all wit there
% maximum at t=12h. The columns are the s-dimensional lieklihood curves.
% This is useful fro getting an idea of what  the likelihood ratio curves
% and Thetas are for a family of curves
% 

figure
hold on
numftimes = size(Mall,1);
for j=1:size(Mall,2)
    M2 = Mall(:,j)';
    M22=squeeze(M2./max(M2));
    tt=24*(1:numftimes)/numftimes;
    N=numftimes;
    Taim=floor(numftimes/2);

    [mx,T]=max(M22);
    Tadj = round(mod(1-Taim+T:N-Taim+T,N))+1;
    plot(tt,M22(Tadj))
    if or(j == 1 , j == 20)
        plot(tt,M22(Tadj),'Linewidth', 2);
    else
        plot(tt,M22(Tadj),'Linewidth', 1)
    end
end
eta=0.35;
epsilon=0.4;
ttt=tt/24;
C=eta*(1+epsilon + cos(2*pi*((ttt-Taim/(N)))));
plot(tt,C,'k','Linewidth', 1)
xlim([0 24])
hold off