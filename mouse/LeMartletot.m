clear all

load('data_hogen')
load('esetfrmaMart_data')

Inx11=[29313 6747 25191 28920 7808 11341 14457 15237 29300 9574 26358];
%% in this dataset the last 2 datapoints should be swapped - they were just read in out of order
organnames={'Adrenal','Aorta','BrownFat','Heart','Kidney','Liver','Lung','SkelMus'};
times = {'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty','twentyone','twentytwo','twentythree','twentyfour'};


clock11Martletot=esetfrmaMartletot_data(Inx11,:);
 
hogen11=esetfrmaHogencore(Inx11,:);

clear clock12Martletotn
for i=1:7;
clock12Martletotn(:,i) = (clock11Martletot(:,i) - mean(clock11Martletot(:,i)))/std(clock11Martletot(:,i));
end

Inx11names={"Arntl",
"Npas2",
"Dbp",
"Per3",
"Nr1d1",
"Per2",
"Nr1d2",
"Tef",
"Ciart",
"Wee1",
"Clock"}

clear clock12hogenstr
for l = 1:24;
    o=0;
          for j=[1 2 4 6 8 9 10 11];
                 o = o+1;
clock12hogenstr.(times{l})(:,o) = hogen11(:,24*(j-1)+l);
          end
    end

clear clock12hogenstrN
for l=1:24;
for y=1:8;
clock12hogenstrN.(times{l})(:,y) = (clock12hogenstr.(times{l})(:,y)- mean(clock12hogenstr.(times{l})(:,y)))/std(clock12hogenstr.(times{l})(:,y));
end
end




clear clock12hogenstrNd test
clock12hogenstrNd=clock12hogenstrN;


test = clock12Martletotn;


clear A U sig V sigma GMModel PROJECTION
 for j = 1:24;
     A.(times{j}) = clock12hogenstrNd.(times{j})*(1/sqrt(8-1)); 
[U.(times{j}),sig.(times{j}),V.(times{j})]=svd(A.(times{j})); 
 end

 
 clear PROJECTION
 for j =1:24;
 for i = 1:24;
PROJECTION.(times{j}).(times{i}) = U.(times{i})'* clock12hogenstrNd.(times{j}); 
PROJECTION.test.(times{i}) = U.(times{i})'*test; 
 end
end
% % 
% cm = [hsv(12); hsv(12)];
%  figure()
%              for i = 1:12;% 1:24 % this is time as well
%                  subplot(3,4,i) 
%             for j =1:24;
%             hold on
%             plot3(PROJECTION.(times{j}).(times{i})(1,:),PROJECTION.(times{j}).(times{i})(2,:),PROJECTION.(times{j}).(times{i})(3,:),'o','Color',cm(j,:),'MarkerFaceColor',cm(j,:),'Markersize',10)
% %             if i==j
% %                plot3(PROJECTION.(times{j}).(times{i})(1,:),PROJECTION.(times{j}).(times{i})(2,:),PROJECTION.(times{j}).(times{i})(3,:),'k.')
% %             end
%             end
%             plot3(PROJECTION.test.(times{i})(1,1:7),PROJECTION.test.(times{i})(2,1:7),PROJECTION.test.(times{i})(3,1:7),'k*','Markersize',10,'Linewidth',2)      
%               xlabel('Latent Variable 1')
%             ylabel('Latent Variable 2')
%             zlabel('Latent Variable 3')
%             grid on
%             title({'Projection with PCs of Time',i})
%          
%              end
%    legend([times,'Original projection', 'Test individual' ])
% % 

clear GMModel
for j =1:24;
 for i = 1:24;%
GMModel.(times{j}).(times{i}) = fitgmdist(PROJECTION.(times{j}).(times{i})(1:3,:)',1);
 end
end




%% Here starts cont
signam = {'one','two','three','four','five','six','seven','eight','nine'};

int = 0.125/2; %interval between data points
intst= 0:int*2:6*4;; %interval between data points

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
%plot3(mus_all(1,:),mus_all(2,:),mus_all(3,:),'k','Linewidth',2)
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
   

clear m2 M2 m4
for r = 1:24;
for q =1:length(Mu_struct.(times{r}));
m2(r,q,:) = mvnpdf(PROJECTION.test.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3));
end
end

m4 = m2;

for r = 1:24;
for p =1:7;
for q=1:769;
 

if m4(r,q,p) < exp(-5);% exp(-5)
    m4(r,q,p) = exp(-5);% exp(-5)
end
    

M2(q,p) =  geomean(m4(:,q,p));
end
end
end

clear L L2 func a a2




eta=0.35;%0.35;
epsil=0.4;
clear dys_new prob prop

M3=sqrt(M2(1:384,:).*M2(385:768,:));
for p = 1:7;
[L2(p), a2(p)] = max(M3(:,p));
end


eta=0.35;
epsil=0.4;
clear dys_new prob prop
for yu = 1:7;
T = find(M3(:,yu)== max(M3(:,yu)));
T=T(1);
if T>192
   Forward = [T+1:384,1:(T-192)]; 
    Backward = [T-192:T-1];
else
     Forward = [T+1:(T+192)];
   Backward = [192+T:384,1:T-1];
end

for c =1:192;
LF(c) = 1/(epsil+1+cos(c/192*pi)).*M3(Forward(c),yu);
end
for c =1:192;
LB(c) =1/(epsil+1+cos((193-c)/192*pi)).*M3(Backward(c),yu);
end

LFN = LF./M3(T,yu);
 LBN = LB./M3(T,yu);

prop = find([LFN LBN]>eta);
prob(yu) = length(prop)/length(M3); 

end



figure(); hold on
plot(linspace(0,24,384),M3)
xlabel('Time')
ylabel('Likelihood')
axis([0 24 0 inf])
box on
grid on

%% this is the order the data is in
actuals = [ 2 6 10 14 18 26 22];

asc = a2/384*24;

asc=mod(asc-6,24);


figure()
hold on
plot(actuals,asc,'r*')
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

figure()
boxplot(prob)
ylabel('Theta')
set(gca,'XTickLabel',{'WT'})

figure()
for i = 1:11;
    subplot(4,3,i)
    hold on
          for j=[1 2 4 6 8 9 10 11 12];       
 plot(18:2:64,hogen11(i,24*(j-1)+1:24*j),'k-');
          end
plot(2+24,clock11Martletot(i,1),'bo') 
plot(6+24,clock11Martletot(i,2),'ro') 
plot(10+24,clock11Martletot(i,3),'mo')
plot(14+24,clock11Martletot(i,4),'go') 
plot(18+24,clock11Martletot(i,5),'b*') 
plot(22+24,clock11Martletot(i,6),'r*') 
plot(26+24,clock11Martletot(i,7),'m*')
axis([18 66 5 12])
title(Inx11names(i))
xlabel('Time')
ylabel('Expression')
end

