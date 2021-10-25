function compare_indiv_loglikelis(inds, lik1,lik2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
L=size(lik1,1);
x=24*(1:L)/L;

hold on
for i = inds
    y1=lik1(:,i);
    y2=lik2(:,i);
    [mx1,T1]=max(y1);
    [mx2,T2]=max(y2);
    y1=y1/mx1;
    y2=y2/mx2;
    plot(x,y1,'b');
    plot(x,y2,'r');
end
hold off
return


