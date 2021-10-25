function  [M M2 m m2] = Lazar
clear all
load('data_hogen')
load('esetfrmaLazar')

Inx11=[29313 6747 25191 28920 7808 11341 14457 15237 29300 9574 26358];
Inx10=[29313 6747 25191 28920 7808 14457 15237 29300 9574 26358];

Inx11=Inx10; % takes out Rev-Erb-Alpha

Inx11names={"Arntl",
"Npas2",
"Per3",
"Dbp",
"Per2",
"Nr1d1",
"Nr1d2",
"Tef",
"Wee1",
"Ciart",
"Clock"};

Inx11names={"Arntl",
"Npas2",
"Per3",
"Dbp",
"Per2",
"Nr1d2",
"Tef",
"Wee1",
"Ciart",
"Clock"}



organnames={'Adrenal','Aorta','BrownFat','Heart','Kidney','Liver','Lung','SkelMus'};
times = {'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty','twentyone','twentytwo','twentythree','twentyfour'};


hogen11=esetfrmaHogencore(Inx11,:);

clock11Lazar=esetfrmaLazar(Inx11,:);


clear clock11hogenstr
for l = 1:24;
    o=0;
          for j=[1 2 4 6 8 9 10 11];
                 o = o+1;
clock11hogenstr.(times{l})(:,o) = hogen11(:,24*(j-1)+l);
          end
    end

clear clock11hogenstrN
for l=1:24;
for y=1:8;
clock11hogenstrN.(times{l})(:,y) = (clock11hogenstr.(times{l})(:,y)- mean(clock11hogenstr.(times{l})(:,y)))/std(clock11hogenstr.(times{l})(:,y));
end
end

clear clock11Lazarn
for i=1:10;
clock11Lazarn(:,i) = (clock11Lazar(:,i) - mean(clock11Lazar(:,i)))/std(clock11Lazar(:,i));
end


clear clock11hogenstrNd test
clock11hogenstrNd=clock11hogenstrN;

for j = 1:24;
test = clock11Lazarn;
end

clear A U sig V sigma GMModel PROJECTION
 for j = 1:24;
     A.(times{j}) = clock11hogenstrNd.(times{j})*(1/sqrt(8-1)); 
[U.(times{j}),sig.(times{j}),V.(times{j})]=svd(A.(times{j})); 
 end

 
 clear PROJECTION
 for j =1:24;
 for i = 1:24;
PROJECTION.(times{j}).(times{i}) = U.(times{i})'* clock11hogenstrNd.(times{j}); 
PROJECTION.test.(times{i}) = U.(times{i})'*test; 
 end
end
% 
cm = [hsv(12); hsv(12)];
 figure()
             for i = 1:12;% 1:24 % this is time as well
                 subplot(3,4,i) 
            for j =1:24;
            hold on
            plot3(PROJECTION.(times{j}).(times{i})(1,:),PROJECTION.(times{j}).(times{i})(2,:),PROJECTION.(times{j}).(times{i})(3,:),'o','Color',cm(j,:),'MarkerFaceColor',cm(j,:),'Markersize',10)
%             if i==j
%                plot3(PROJECTION.(times{j}).(times{i})(1,:),PROJECTION.(times{j}).(times{i})(2,:),PROJECTION.(times{j}).(times{i})(3,:),'k.')
%             end
            end
            plot3(PROJECTION.test.(times{i})(1,1:5),PROJECTION.test.(times{i})(2,1:5),PROJECTION.test.(times{i})(3,1:5),'k*','Markersize',10,'Linewidth',2)      
              plot3(PROJECTION.test.(times{i})(1,6:10),PROJECTION.test.(times{i})(2,6:10),PROJECTION.test.(times{i})(3,6:10),'ko','Markersize',10,'Linewidth',2)
              xlabel('Latent Variable 1')
            ylabel('Latent Variable 2')
            zlabel('Latent Variable 3')
            grid on
            title({'Projection with PCs of Time',i})
         
             end
   legend([times,'Original projection', 'Test individual' ])
% 

clear GMModel
for j =1:24;
 for i = 1:24;%
GMModel.(times{j}).(times{i}) = fitgmdist(PROJECTION.(times{j}).(times{i})(1:3,:)',1);
 end
end




%% Here starts cont
signam = {'one','two','three','four','five','six','seven','eight','nine'};

int = 0.125/2; %interval between data points
%intst= 0:int*2:6*4;; %interval between data points

clear Mu_struct Sigma_struct
for var = 1:24;
clear mus sigmas pp ppm mus_all Sigma_mat
for j=1:24;
mus(:,j) = GMModel.(times{j}).(times{var}).mu;
sigmas(:,j) =   GMModel.(times{j}).(times{var}).Sigma(:);
end
mus(:,25) = GMModel.(times{1}).(times{var}).mu;
sigmas(:,25) =   GMModel.(times{1}).(times{var}).Sigma(:);
for i = 1:9;
pp1.(signam{i}) =pchip(0:2:48,[sigmas(i,:)']);
end
%% matches first and second derivatives at end points
ppm = csape([0:2:48]',[mus(1,:)',mus(2,:)',mus(3,:)']','periodic');
mus_all =ppval(ppm,0:int:48); % why is 24 not enough? That's how the function was built!!
% 
% hold on
% subplot(2,3,var)
plot3(mus_all(1,:),mus_all(2,:),mus_all(3,:),'k','Linewidth',2)
% grid on

for i = 1:9;
Sigma_mat(i,:) = ppval(pp1.(signam{i}),0:int:48); %
end

for q = 1:length(mus_all);
if find(eig(reshape(Sigma_mat(:,q),3,3)) <=0)
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
   

clear m m2 M M2
for r = 1:24;
for q =1:length(Mu_struct.(times{r}));
m(r,q,:) = log(mvnpdf(PROJECTION.test.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3)));  
m2(r,q,:) = mvnpdf(PROJECTION.test.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3));
end
end

m3=m;  
m4 = m2;

for r = 1:24;
for p =1:10;
for q=1:769;
 
if m3(r,q,p) < -5;% -5;
    m3(r,q,p) = -5; %-5
end

if m4(r,q,p) < exp(-5);% exp(-5)
    m4(r,q,p) = exp(-5);% exp(-5)
end
    
M(q,p) = mean(m3(:,q,p));
M2(q,p) =  geomean(m4(:,q,p));
end
end
end

clear L L2 func a a2
for p = 1:10;
[L(p), a(p)] = max(M(:,p));
end

for p = 1:10;
[L2(p), a2(p)] = max(M2(:,p));
end
% 
% for p = 1:10;
% func(p) = max(M2(:,p))./trapz(M2(:,p))*100;
% end



figure(); hold on
plot(18:0.05985:64,M2(:,1:5),'b')
plot(18:0.05985:64,M2(:,6:10),'r')
xlabel('Time')
ylabel('Likelihood')
%axis([18 64 0 1.0])
box on
grid on


a/769*46 +18




% figure()
% boxplot([func(1:5);func(6:10)]')
% ylabel('Confidence Measure')
% set(gca,'XTickLabel',{'WT','KO'})





figure()
for i = 1:10;
    subplot(4,3,i)
    hold on
          for j=[1 2 4 6 8 9 10 11 12];       
 plot(18:2:64,hogen11(i,24*(j-1)+1:24*j),'b-');
          end
plot(34,clock11Lazar(i,1:5),'k*') 
plot(34,clock11Lazar(i,6:10),'ro') 
axis([18 64 -inf inf])
title(Inx11names(i))
xlabel('Time')
ylabel('Expression')
end


eta=0.35;%0.35;
epsil=0.4;
clear dys_new prob prop


clear M3 a3 L3
M3=sqrt(M2(1:384,:).*M2(385:768,:));
for p = 1:10;
[L3(p), a3(p)] = max(M3(:,p));
end

clear prob LF LB LFN LBN prop
for yu = 1:10;
T = find(M3(:,yu)== max(M3(:,yu)));
T=T(1);
if T>192%T>96
    Forward = [T+1:384,1:(T-192)];
    %Forward = [T+1:193,1:(T-97)];
    Backward = [T-192:T-1];
else
     Forward = [T+1:(T+192)];
    Backward = [193+T:384,1:T-1];
end
%% check
% figure(); hold on;  plot(Forward,M2(Forward,yu));
%  hold on; plot(Backward,M2(Backward,yu)) 
for c =1:192;
LF(c) = 1/(epsil+1+cos(c/192*pi)).*M3(Forward(c),yu);
end
for c =1:192;
LB(c) =1/(epsil+1+cos((193-c)/192*pi)).*M3(Backward(c),yu);
end
% LFN = M2(T,yu)./LF;
 %LBN = M2(T,yu)./LB;
LFN = LF./M3(T,yu);
 LBN = LB./M3(T,yu);

prop = find([LFN LBN]>eta);
prob(yu) = length(prop)/length(M3); 
%dys_new(yu) =log10(sum([LFN LBN]));
end
prob;


figure(); hold on
plot(linspace(18,18+24,length(M3)),M3(:,1:5),'b')
plot(linspace(18,18+24,length(M3)),M3(:,6:10),'r')
xlabel('Time')
ylabel('Likelihood')
%axis([18 64 0 1.0])
box on
grid on

a3*0.0627 +18




figure()
hold on
boxplot(prob(:),[ones(1,5) 2*ones(1,5)])
ylabel('\Theta')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'XTickLabel',{'WT', 'KO'})
grid on

