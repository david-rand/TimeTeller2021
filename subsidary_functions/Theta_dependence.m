function Thetas = Theta_dependence(likelis, etamin, etamax)

epsil=0.4;
etas=etamin:(etamax-etamin)/10:etamax;
Thetas=zeros(length(etas),size(likelis,2));
for j=1:size(likelis,2)
    for i=1:length(etas)
        likeli=likelis(:,j);
        Thetas(i,j) = Theta_from_likelihood( likeli,etas(i),epsil );
    end
end
% plot(etas,Thetas,'o');

function Theta = Theta_from_likelihood( likeli,eta,epsil )
% Extarcts Theta values from single likelihood functions
%   Detailed explanation goes here

num_interp_times=length(likeli);
clear dys_new prob prop Forward Backward LFN LBN LF LB
first_half = floor(num_interp_times/2);
    %figure(); plot(M2(:,yu))
T = find(likeli== max(likeli));
T=T(1);
if T>first_half
    Forward = [T+1:num_interp_times,1:(T-(first_half+1))];
    Backward = [T-first_half:T-1];
else
    Forward = [T+1:(T+first_half)];
    Backward = [(first_half+1)+T:num_interp_times,1:T-1];
end

for c =1:first_half;
    LF(c) = 1/(epsil+1+cos(c/96*pi)).*likeli(Forward(c));
end
for c =1:first_half;
    LB(c) =1/(epsil+1+cos(((first_half+1)-c)/first_half*pi)).*likeli(Backward(c));
end

LFN = LF./likeli(T);
LBN = LB./likeli(T);

%figure(2); plot(Backward,LBN,'*'); hold on; plot(Forward,LFN,'*')
prop = find([LFN LBN]>eta);
Theta = length(prop)/num_interp_times; 
% gives Theta
return

