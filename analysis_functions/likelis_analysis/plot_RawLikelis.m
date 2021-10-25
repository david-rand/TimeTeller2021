function plot_RawLikelis(RawLikelis, Likelis,ind,logthresh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x=(1:193)*24/193;
Z=max(RawLikelis(:,:,ind)',logthresh);
Y=sum(Z')/6;
plot(x,Z,'LineWidth',2)
%Y=Y - max(Y);
hold on
plot(x,Y,'k','LineWidth',2)
hold off
end

