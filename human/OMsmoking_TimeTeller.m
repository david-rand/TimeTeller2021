function OM_smoking_what_time_july17

clear all
load('esetfrmaSMOM_data')
load('frmaOM')
load('Inx16') %Indexes of the training genes
times = {'am8','midday','pm4','pm8','midnight','am4'};
Inx = Inx16;
Inx([15])=[]; %the 15th probe is unusable here - expression is too high
Inx([11])=[]; %t
clear SMOMclock clock_OM
SMOMclock = esetfrmaSMOM_data(Inx,:); %new
clock_OM = frmaOM(Inx,:);

h=length(Inx);

figure()
hold on
for i = 1:h;
subplot(4,4,i)
% set(p, 'position',[0 0 0.4 0.4]);
hold on
boxplot(SMOMclock(i,:),[ones(1,40),2*ones(1,39)],'boxstyle','filled') % control
plot(repmat(3:8,1,10),clock_OM(i,:),'k*')
h2=area([0.5 8.5],[max(clock_OM(i,:)),max(clock_OM(i,:))],'BaseValue',[min(clock_OM(i,:))])
h2.FaceColor = [0.25 0.25 0.25];
h2.FaceAlpha = 0.1;
 axis([0.5 8.5 3 11.5])
%title(clock16OMnames{i})
 xticks([ 1 2 3 4 5 6 7 8 ])
xticklabels({'Non Smoker','Smoker', ' ','8am','12d','4pm','8pm','12n','4am'})
%set(gca,'FontSize',10)
%xtickangle(90)
 %axis equal;
%   xlim([0.5, 16.5]);
%  ylim([3, 11.5]);
end

clear clock_OMN SMOMclockN

for n=1:79;
    for m = 1:h;
        M=mean(SMOMclock(:,n));
        S = std(SMOMclock(:,n));
SMOMclockN(m,n) = (SMOMclock(m,n)-M )/(S);
    end
end

for n=1:60;
    for m = 1:h;
      M=mean(clock_OM(:,n));
        S = std(clock_OM(:,n));
 clock_OMN(m,n) = (clock_OM(m,n)-M)/(S);
    end
end

figure()
hold on
for i = 1:h;
subplot(4,4,i)
% set(p, 'position',[0 0 0.4 0.4]);
hold on
boxplot(SMOMclockN(i,:),[ones(1,40),2*ones(1,39)],'boxstyle','filled') % control
plot(repmat(3:8,1,10),clock_OMN(i,:),'k*')
h2=area([0.5 8.5],[max(clock_OMN(i,:)),max(clock_OMN(i,:))],'BaseValue',[min(clock_OMN(i,:))])
h2.FaceColor = [0.25 0.25 0.25];
h2.FaceAlpha = 0.1;
 axis([0.5 8.5 -inf inf])
%title(probenames_16{Inx(i)})
 xticks([ 1 2 3 4 5 6 7 8 ])
xticklabels({'Non Smoker','Smoker', ' ','8am','12d','4pm','8pm','12n','4am'})
%set(gca,'FontSize',10)
%xtickangle(90)
 %axis equal;
%   xlim([0.5, 16.5]);
%  ylim([3, 11.5]);
end



clear clockOMstr
for j = 1:6;
clockOMstr.(times{j}) = clock_OMN(:,j:6:60);
end

%% have to clear old structures to update size
 clear A U sig V sigma

 %% Do the SVD
 for j = 1:6;
 A.(times{j}) =clockOMstr.(times{j})*(1/sqrt(6-1)); 
[U.(times{j}),sig.(times{j}),V.(times{j})]=svd(A.(times{j})); 
sigma.(times{j})=diag(sig.(times{j})); %singular values already in order
sigma_norm.(times{j}) = sigma.(times{j})/sum(sigma.(times{j}));
 end
 
 %% minimum variance explained by 2PCs is ~82% (better than 3 components in mouse model) third component messes up prediction

 clear test
 test =  SMOMclockN;
 clear PROJECTION
 for j =1:6;
 for i = 1:6;
PROJECTION.(times{j}).(times{i}) = U.(times{i})'*clockOMstr.(times{j}); 
PROJECTION.test.(times{i}) = U.(times{i})'*test; 
 end
 end
 
 mpp2 = hsv(6);
 figure()
             for i = 1:6;% 1:24 % this is time as well
                 subplot(2,3,i)
            for j =1:6;%[1 3 6 9 12 15 18 21 24];
            hold on
            plot3(PROJECTION.(times{j}).(times{i})(1,:),PROJECTION.(times{j}).(times{i})(2,:),PROJECTION.(times{j}).(times{i})(3,:),'o','Color',mpp2(j,:),'MarkerFaceColor',mpp2(j,:),'Markersize',10)
            end
            plot3(PROJECTION.test.(times{i})(1,1:40),PROJECTION.test.(times{i})(2,1:40),PROJECTION.test.(times{i})(3,1:40),'k+','Markersize',10)
            plot3(PROJECTION.test.(times{i})(1,41:79),PROJECTION.test.(times{i})(2,41:79),PROJECTION.test.(times{i})(3,41:79),'ko','Markersize',10)
              xlabel('PC1')
            ylabel('PC2')
            zlabel('PC3')
            grid on
            title(['' ])
         
             end
   legend([times, 'Healthy','Smoker' ])
   
   clear GMModel M M2 m m2 
for j =1:6;
 for i = 1:6;%
GMModel.(times{j}).(times{i}) = fitgmdist(PROJECTION.(times{j}).(times{i})(1:3,:)',1);
 end
end


%% Here starts cont
signam = {'one','two','three','four','five','six','seven','eight','nine', 'ten'};

int = 0.125/2; %interval between data points
intst= 0:int*2:6*4;; %interval between data points

clear Mu_struct Sigma_struct M M2 m m2 m3 m4
for var = 1:6;
clear mus sigmas pp ppm mus_all Sigma_mat 
for j=1:6;
mus(:,j) = GMModel.(times{j}).(times{var}).mu;
sigmas(:,j) =   GMModel.(times{j}).(times{var}).Sigma(:);
end
mus(:,7) = GMModel.(times{1}).(times{var}).mu;
sigmas(:,7) =   GMModel.(times{1}).(times{var}).Sigma(:);
for i = 1:9;
pp1.(signam{i}) =pchip(0:4:24,[sigmas(i,:)']);
end
%% matches first and second derivatives at end points
ppm = csape([0:4:24]',[mus(1,:)',mus(2,:)',mus(3,:)']','periodic');
mus_all =ppval(ppm,0:2*int:24); % why is 24 not enough? That's how the function was built!!

hold on
subplot(2,3,var)
plot3(mus_all(1,:),mus_all(2,:),mus_all(3,:),'k','Linewidth',2)
grid on

for i = 1:9;
Sigma_mat(i,:) = ppval(pp1.(signam{i}),0:int*2:24); %
end

for q = 1:length(mus_all);
if find(eig(reshape(Sigma_mat(:,q),3,3)) <=0.0001)
A=reshape(Sigma_mat(:,q),3,3);
B = (A + A')/2; 
[U5,Sigma5] = eig(B); 
Ahat5 = U5*max(Sigma5,0)*U5';
Ahat5=Ahat5+0.001*eye(3);
Sigma_mat(:,q) =Ahat5(:);
end

Sigma_struct.(times{var})(:,q) = Sigma_mat(:,q);
Mu_struct.(times{var})(:,q) = mus_all(:,q);
end

end
   
%% Likelihoods
for r = 1:6;
for q =1:length(Mu_struct.(times{r}));
m2(r,q,:) = mvnpdf(PROJECTION.test.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3));
end
end

m4 = m2;

for r = 1:6;
for p =1:79;
for q=1:193;
 
if m4(r,q,p) < exp(-5);% exp(-5)
    m4(r,q,p) = exp(-5);% exp(-5)
end
    
M2(q,p) =  geomean(m4(:,q,p));

end
end
end

clear L2 a2
for p = 1:79;
[L2(p), a2(p)] = max(M2(:,p));
end

a2 = a2/193 *24;



for r =1:6;
lengths(r) = length(Mu_struct.(times{r}));
end

for q=1:min(lengths);
M2(q) =  geomean(m4(:,q));
end


h1 =hist(a2(1:40),[0:1:23])
h2 =hist(a2(41:79),[0:1:23])
figure()
hold on
grid on
H=bar([0:1:23],[h1.' h2.'])
axis([0 24 0 inf])
xticks([0:1:23])
xtickangle([90])
set(gca,'XTickLabel',{'8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00','00:00','01:00','02:00','03:00','04:00','05:00','06:00','7:00'})
H(1).Parent.Parent.Colormap = [0 0 0; [1 0 0]];
xlabel('Maximum Likelihood Estimate')
ylabel('# Estimations')
legend('NonSmoker','Smoker')





figure()
hold on
plot(a2(1:40),L2(1:40),'k+','Markersize',12,'Linewidth',3)
plot(a2(41:79),L2(41:79),'rx','Markersize',12,'Linewidth',3)
hold on
for p = 1:40;
hold on
%plot(m2')
a3 = plot(intst,M2(:,p),'k')
a3.Color(4) = 0.3;
end
for p = 41:79;
hold on
%plot(m2')
a4=plot(intst,M2(:,p),'Color',[0.6 0 0])
a4.Color(4) = 0.4;
end
set(gca,'XTick' ,0:4:20)
set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00'})
legend('NonSmoker','Smoker')




figure()
hold on
plot(a2(1:40),L2(1:40),'k+','Markersize',12,'Linewidth',3)
plot(a2(41:79),L2(41:79),'rx','Markersize',12,'Linewidth',3)
for p = 1:40;
hold on
%plot(m2')
a3 = plot(intst,M2(:,p),'k')
a3.Color(4) = 0.3;
end
for p = 41:79;
hold on
%plot(m2')
a4=plot(intst,M2(:,p),'Color',[0.6 0 0])
a4.Color(4) = 0.4;
end
hold on
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00','08:00'})
legend('NonSmoker','Smoker')
xlabel('Time')
ylabel('Log Likelihood')
box on
grid on
axis([0 24 0 0.85])


clear dys_new prop prob
epsil=0.4; 
eta = 0.35;
for yu = 1:79; 

clear LBN LFN Backward Forward LB LF
T = find(M2(:,yu)== max(M2(:,yu)));
T=T(1);
if T>96
    Forward = [T+1:193,1:(T-97)];
    Backward = [T-96:T-1];
else
     Forward = [T+1:(T+96)];
    Backward = [97+T:193,1:T-1];
end
for c =1:96;
LF(c) = 1/(epsil+1+cos(c/96*pi)).*M2(Forward(c),yu);
end
for c =1:96;
LB(c) =1/(epsil+1+cos((97-c)/96*pi)).*M2(Backward(c),yu);
end


LFN =LF./M2(T,yu);
LBN =LB./M2(T,yu);

%figure(2); plot(Backward,LBN,'-'); hold on; plot(Forward,LFN,'-')

prop = find([LFN LBN]>eta);
prob(yu) = length(prop)/193; 
end

%prob







figure()
hold on
plot(a2(1:40),prob(1:40),'k+','Markersize',12,'Linewidth',3)
plot(a2(41:79),prob(41:79),'rx','Markersize',12,'Linewidth',3)
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00','08:00'})
legend('NonSmoker','Smoker')
xlabel('Time')
ylabel('Log Likelihood')
box on
grid on
axis([0 24 0 0.85])


figure()
boxplot([prob(1:40) ,prob(41:79)],[ones(1,40), 2*ones(1,39)])%,dyssmoker])
set(gca,'XTickLabel',{'Normal','Smoker'})







figure()
hold on
boxplot([prob(1:40) ,prob(41:79)],[ones(1,40), 2*ones(1,39)])%,dyssmoker])
plotSpread([prob'],'distributionIdx',[ones(1,40), 2*ones(1,39)],'distributionMarkers',{'.','.'},'distributionColors',{'k','k'})
set(gca,'XTickLabel',{'Non-Smoker n=40','Smoker n=39'})
set(findobj(gca,'type','line'),'linew',2)


figure()
hold on
plotSpread([prob'],'distributionIdx',[ones(1,40), 2*ones(1,39)],'distributionMarkers',{'.','.'},'distributionColors',{'k','k'})
boxplot([prob(1:40) ,prob(41:79)],[ones(1,40), 2*ones(1,39)])%,dyssmoker])
set(gca,'XTickLabel',{'Non-Smoker n=40','Smoker n=39'})
pval_lohav =ranksum(prob(1:40),prob(41:79))
title({ sprintf('Significance p=%0.4f' ,pval_lohav)})
box on
set(findobj(gca,'type','line'),'linew',2)
grid on
ylabel('\Theta')

hold on
meanMhealthy = mean(M(:,1:40)');
meanMsmoker = mean(M(:,41:79)');

meanM2healthy = mean(M2(:,1:40)');
meanM2smoker = mean(M2(:,41:79)');

figure()
hold on
plot(meanM2healthy,'k')
plot(meanM2smoker,'color',[0.6 0 0])



figure()
hold on
subplot(2,2,1)
for p = 1:20;
    hold on
plot(linspace(0,24,193),M2(:,p),'k')
end
title('Non-smoker. Female')
%set(gca,'XTick',1:6)
%set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00'})
subplot(2,2,2)
for p = 21:40;
hold on
plot(linspace(0,24,193),M2(:,p),'k')
end
title('Non-smoker. Male')
%set(gca,'XTick',1:6)
%set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00'})
subplot(2,2,3)
for p = 41:59;
hold on
plot(M2(:,p),'Color',[0.6 0 0])
end
set(gca,'XTick',1:6)
set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00'})
title('Smoker. Female')
subplot(2,2,4)
for p = 60:79;
hold on
plot(M2(:,p),'Color',[0.6 0 0])
end
title('Smoker. Male')
set(gca,'XTick',1:6)
set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00'})
suptitle('Likelihoods')





