function sriuk_TimeTeller

clear all
load('esetfrmasriukOM_data');
load('Inx16');
load('frmaOM');
Inx = Inx16;
Inx(15)=[];
sriukclock = esetfrmasriukOM_data(Inx,:);
esetfrmaOM= frmaOM(Inx,:);
h=length(Inx);
times = {'am8','midday','pm4','pm8','midnight','am4'};




for n=1:8;
    for m = 1:h;
       s=std(sriukclock(:,n));
        me= mean(sriukclock(:,n));
  sriukclockN(m,n) = (sriukclock(m,n)-me )/(s);
    end
end
for n=1:60;
    for m = 1:h;
        s=std(esetfrmaOM(:,n));
        me= mean(esetfrmaOM(:,n));
esetfrmaOM_dataN(m,n) = (esetfrmaOM(m,n)-me )/(s);
    end
end





for j = 1:6;
esetfrmaOM_datastr.(times{j}) = esetfrmaOM_dataN(:,j:6:60);
end


%% have to clear old structures to update size
 clear AN UN sigN VN sigmaN sigma_normN

 %% Do the SVD
 for j = 1:6;
 AN.(times{j}) =esetfrmaOM_datastr.(times{j})*(1/sqrt(6-1)); % should this be 8 genes or 9 organs? hmmmm (doens't really matter)
[UN.(times{j}),sigN.(times{j}),VN.(times{j})]=svd(AN.(times{j})); 
sigmaN.(times{j})=diag(sigN.(times{j})); %singular values already in order
sigma_normN.(times{j}) = sigmaN.(times{j})/sum(sigmaN.(times{j}));
 end
 
 %% minimum variance explained by 2PCs is ~82% (better than 3 components in mouse model) third component messes up prediction
 figure()
 for j = 1:6;
 plot(cumsum(sigma_normN.(times{j})),'k*-')
 hold on
 end
 

clear PROJECTIONNN
 for j =1:6;
 for i = 1:6;
PROJECTIONN.(times{j}).(times{i}) = UN.(times{i})'*esetfrmaOM_datastr.(times{j}); 
PROJECTIONN.test.(times{i}) = UN.(times{i})'* sriukclockN; 
 end
 end
 mpp2=hsv(6);
              figure()
             for i = 1:6;% 1:24 % this is time as well
                 subplot(2,3,i)
            for j =1:6;%[1 3 6 9 12 15 18 21 24];
            hold on
            plot3(PROJECTIONN.(times{j}).(times{i})(1,:),PROJECTIONN.(times{j}).(times{i})(2,:),PROJECTIONN.(times{j}).(times{i})(3,:),'o','Color',mpp2(j,:),'MarkerFaceColor',mpp2(j,:),'Markersize',10)
            end
            plot3(PROJECTIONN.test.(times{i})(1,1:5),PROJECTIONN.test.(times{i})(2,1:5),PROJECTIONN.test.(times{i})(3,1:5),'k*','Markersize',10)
             plot3(PROJECTIONN.test.(times{i})(1,6:8),PROJECTIONN.test.(times{i})(2,6:8),PROJECTIONN.test.(times{i})(3,6:8),'ko','Markersize',10)
            %  plot3(Mu_struct.(times{i})(1,:),  Mu_struct.(times{i})(2,:), Mu_struct.(times{i})(3,:), 'k')
            xlabel('PC1')
            ylabel('PC2')
            zlabel('PC3')
            grid on
            axis equal
            title(['' ])
         
             end


              
   

clear GMModelN M M2 m m2 
for j =1:6;
 for i = 1:6;%
GMModelN.(times{j}).(times{i}) = fitgmdist(PROJECTIONN.(times{j}).(times{i})(1:3,:)',1);
 end
end


%% Here starts cont
signam = {'one','two','three','four','five','six','seven','eight','nine', 'ten'};

int = 0.125/2; %interval between data points
intst= 0:int*2:6*4;

clear Mu_struct Sigma_struct

for var = 1:6;
clear mus sigmas pp ppm mus_all Sigma_mat 
for j=1:6;
mus(:,j) = GMModelN.(times{j}).(times{var}).mu;
sigmas(:,j) =   GMModelN.(times{j}).(times{var}).Sigma(:);
end
mus(:,7) = GMModelN.(times{1}).(times{var}).mu;
sigmas(:,7) =   GMModelN.(times{1}).(times{var}).Sigma(:);

for i = 1:9;
pp1.(signam{i}) =pchip(0:4:24,[sigmas(i,:)']);
end

%% matches first and second derivatives at end points
ppm = csape([0:4:24]',[mus(1,:)',mus(2,:)',mus(3,:)']','periodic');

%ppm = pchip(0:4:24,[mus(1,:)',mus(2,:)',mus(3,:)']');

% 
% figure()
% hold on
% for i = 1:9;
% fnplt(pp1.(signam{i}),'b')%'Linewidth',2)
% %fnplt(pp.(signam{i}),'k-.')
% end
% hold on
% plot(0:4:24,sigmas(1,:),'rx','Markersize',8,'Linewidth',3)
% legend('pchip','spline','Data Fits')


mus_all =ppval(ppm,0:2*int:24); % why is 24 not enough? That's how the function was built!!

%mus_all(:,find(mus_all(1,:)>25))=[];
% 

% hold on
% subplot(2,3,var)
% plot3(mus_all(1,:),mus_all(2,:),mus_all(3,:),'k','Linewidth',2)
% % plot3(mus(1,:),mus(2,:),mus(3,:),'ro','Markerface','r')
% xlabel('LatentVar1')
% ylabel('LatentVar2')
% zlabel('LatentVar3')
% hold on
% grid on
% legend('csape fit','Latent Variable Means')
% 

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


%%%
clear m m2 m3 m4 M M2 M3 M4 M6 M5 lengths

for r = 1:6;
for q =1:length(Mu_struct.(times{r}));
m2(r,q,:) = mvnpdf(PROJECTIONN.test.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3));
end
end

m4 = m2;

for r = 1:6;
for p =1:8;
for q=1:193;

if m4(r,q,p) < exp(-10);% exp(-5)
    m4(r,q,p) = exp(-10);% exp(-5)
end

M2(q,p) =  geomean(m4(:,q,p));

end
end
end


clear a a2 L2 
for p = 1:8;
[L2(p), a2(p)] = max(M2(:,p));
end


figure()
hold on
plot(linspace(0,24,193),M2(:,1:5),'b','Linewidth',2)
plot(linspace(0,24,193),M2(:,6:8),'g','Linewidth',2)
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00','8:00'})
ylabel('Likelihood')
xlabel('Time')
box on
grid on
legend({'UK','Sri'})




a = a2/193 *24;


eta=0.35;
epsil=0.4;


clear dys_new prop prob

for yu = 1:8; %times?
%figure(); plot(M2(:,yu))
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
%dys_new(yu) =log10(sum([LFN LBN]));
end

prob

figure()
hold on
plot(a(1:5),prob(1:5),'b*','Linewidth',2)
plot(a(6:8),prob(6:8),'g*','Linewidth',2)
axis([0 24 0 0.3])
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'8:00','12:00','16:00','20:00','00:00','04:00','8:00'})
ylabel('\Theta')
xlabel('Time')
box on
grid on