function  mouse_training_TimeTeller

% This function uses a leave on organ out approach to contruct a
% probability model using the data on the other organs and the using this
% model to calculate the likelihood curves for the left out organ. It then
% plots these curves for all organs and plots the estimated time for each
% organ-time data against its true time.

% History
% 2016 Used for Denise Thesis
% 2018-9 Used for paper.
% Last check: May 2019
% Further commenting: May 2019

% this is the code for leave one organ out validation
load('zhang11_12organs_data')
% a 11x288 double array, 11 tissues and 

organnames={'Adrenal','Aorta','BrownFat','Heart','Kidney','Liver','Lung','SkelMus'};
times = {'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty','twentyone','twentytwo','twentythree','twentyfour'};

clear clock12hogenstr
for tms = 1:24; % times
    counter=0;
    for j=[1 2 4 6 8 9 10 11]; % 8 tissues
        counter = counter+1; % allows for possibility that some js are missing
        clock11hogenstr.(times{tms})(:,counter) = zhang11_12organs(:,24*(j-1)+tms);
        % clock11hogenstr.one is a 11x8 arrat, 11 probe expressions times 8 timepoints
    end
end

clear clock12hogenstrN
% gene vectors are normalised
for tms=1:24;
    for tis=1:8;
        clock11hogenstrN.(times{tms})(:,tis) = (clock11hogenstr.(times{tms})(:,tis)- mean(clock11hogenstr.(times{tms})(:,tis)))/std(clock11hogenstr.(times{tms})(:,tis));
    end
end


% ===================================================================
% start of organ loop
clear a a2 L L2 M M2 func Mall
for organ =1:8;
    clear clock11hogenstrNd test
    clock11hogenstrNd=clock11hogenstrN;

    for j = 1:24;
        test(j,:) = clock11hogenstrNd.(times{j})(:,organ);
        % test is 24x11 double times vs genes for organ
    end
    clear A U sig V sigma GMModel PROJECTION
    i=organ;
    for j = 1:24;
        if  i == organ; 
            clock11hogenstrNd.(times{j})(:,i)=[]; 
        end
    % take the organ data out from the probes panel for leave-one-out
    % validation
    A.(times{j}) = clock11hogenstrNd.(times{j});%*(1/sqrt(8-1)); 
    [U.(times{j}),sig.(times{j}),V.(times{j})]=svd(A.(times{j})); 
    % U.one is the projection at time one and sigma and V are the
    % corresponding outputs from the SVD, similarly for other times
    sigma.(times{j})=diag(sig.(times{j})); %singular values already in order
    sigma_norm.(times{j}) = sigma.(times{j})/sum(sigma.(times{j}));
    % these are the singular values and the fractional power/variance
    end

    % PROJECTION.(times{j}).(times{i})(1:3,:)' is projection into 3d
    % The fields of PROJECTION are the times and, for example, PROJECTION.one has the
    % same fields. Then e.g. PROJECTION.one.two is a 11×8 double obtained by
    % using the am8 projection to project the two data into R^3.
    clear PROJECTION
    for j =1:24;
        for i = 1:24;
            PROJECTION.(times{j}).(times{i}) = U.(times{i})'* clock11hogenstrNd.(times{j});
            % this is the data for time t_j projectred using the projection
            % from time t_i
            PROJECTION.test.(times{i}) = U.(times{i})'*test'; 
            % this is the projection by the t_i projector of tissue
            % data test(j,:) i.e. the probe vector at the jth time point of
            % organ
        end
    end


% GMModel = fitgmdist(X,k,Name,Value) returns a Gaussian mixture distribution model 
% with additional options specified by one or more Name,Value pair arguments.
% Then e.g. GMModel.one.two is a fit of a multivariate normal to
% PROJECTION.one.two in 3d
% PROJECTION..one.two uses the two projection to project the one data into R^3.
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
            % mus(:,j) is the mean of the t_j data points after projection by
            % the t_var projector
            sigmas(:,j) =   GMModel.(times{j}).(times{var}).Sigma(:);
        end
        mus(:,25) = GMModel.(times{1}).(times{var}).mu;
        sigmas(:,25) =   GMModel.(times{1}).(times{var}).Sigma(:);
        % as above but covariance matrices
        for i = 1:9;
            pp1.(signam{i}) =pchip(0:2:48,[sigmas(i,:)']);
        end
        % pp = pchip(x,y) returns a piecewise polynomial structure for use with ppval and the spline utility unmkpp.
        % p=pchip(x,y,xq) returns a vector of interpolated values p 
        % corresponding to the query points in xq. The values of p 
        % are determined by shape-preserving piecewise cubic interpolation of x and y.
    %% matches first and second derivatives at end points
        ppm = csape([0:2:48]',[mus(1,:)',mus(2,:)',mus(3,:)']','periodic');
        mus_all =ppval(ppm,0:int:48); % why is 24 not enough? That's how the function was built!!
        % pp = pchip(x,y) returns a piecewise polynomial structure for use with ppval and the spline utility unmkpp.
        % p=pchip(x,y,xq) returns a vector of interpolated values p 
        % corresponding to the query points in xq. The values of p 
        % are determined by shape-preserving piecewise cubic interpolation of x and y.

        for i = 1:9;
            Sigma_mat(i,:) = ppval(pp1.(signam{i}),0:int:48);
            % this is a 9 x N array each of the N entries a 3x3 matrix. It
            % extends the covariance matrices to all time
        end

        % checking for non positive eigenvalues (used rarely, but necessary)
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
    % end of var loop
    % ===================================================================
   

    clear m m2
    for r = 1:24;
        for q =1:length(Mu_struct.(times{r}));
            %m(r,q,:) = log(mvnpdf(PROJECTION.test.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3)));  
            m2(r,q,:) = mvnpdf(PROJECTION.test.(times{r})(1:3,:)',Mu_struct.(times{r})(:,q)',reshape(Sigma_struct.(times{r})(:,q) ,3,3));
            % y = mvnpdf(X,MU,SIGMA) returns the density of the multivariate 
            % normal distribution with mean MU and covariance SIGMA, evaluated at each row of X.
            % this gives the likelihood of the test vectors aling the
            % interpolated times i.e. the likelihood functions for all 6 of
            % the individuals probe vectors
        end
    end
    % This gives curly likelihood curve L_\ell(t|test)

    for p =1:24;
        for q=1:769;
            if m2(r,q,p) < exp(-5);% exp(-5)
                m2(r,q,p) = exp(-5);% exp(-5)
            end
        M2(organ,q,p) =  geomean(m2(:,q,p));
        % M2(organ,:,p) is the lilelihood curve for organ at the pth
        % timepoint
        end
    end
end
% end of organ loop
% ===================================================================

actuals = linspace(0,46,24)';

clear a2 a3 a4 a5 L2 L3 L4 L5
for organ=1:8;
    for p = 1:24;
    [L(p,organ), a(p,organ)] = max(reshape(M2(organ,:,p),769,1)); % PDF multiple
    % value at maximum likelihood and its index
    end
end

a=a/769*48; %converts indices into times

% following figure plots all the likelihood curves at all 12 pairs of times
figure()
hold on
% plot likelihood curves
for organ=1:8;
    for j =1:12;
        hold on
        subplot(4,3,j)
        %p=area(linspace(0,46,769),reshape(M2(organ,:,j),769,1));
        plot(linspace(0,46,769),reshape(M2(organ,:,j),769,1),'b');
    end
    for j =1:12;
        hold on
        subplot(4,3,j)
        plot(linspace(0,46,769),reshape(M2(organ,:,j+12),769,1),'b');
        %q=area(linspace(0,46,769),reshape(M2(organ,:,j+12),769,1));
        %alpha(q,0.2)
    end
end
% plot vertical red lines and labels
for j =1:12;
    hold on
    subplot(4,3,j)
    plot([(j-1)*2,(j-1)*2],[0,6],'r-','Linewidth',3)
    plot([(j+11)*2,(j+11)*2],[0,6],'r-','Linewidth',3)
    axis([0 46 0 max(max(max(M2(:,:,[j j+12]))))])
    set(gca,'xtick',[0:4:48])
    set(gca,'XTickLabel',{'18','22','26','30','34','38','42','46','50','54','58','62','66'})
    xlabel('CT')
    ylabel('Likelihood')
end


% following figure plots estimated time T against actual time
figure()
hold on
for t = 1:24;
    for organ= 1:8;
        hold on
        plot(actuals(t),a(t,organ),'b*','Markersize',6,'Linewidth',1)
    end
end
plot([-1 49], [-1 49],'k-')
plot([-1 25], [23 49],'k-')
plot([23 49], [-1 25],'k-')
plot([-1 1], [47 49],'k-')
plot( [47 49],[-1 1],'k-')

grid on
box on
axis([ -1 48 -1 48])
xlabel('Real CT')
set(gca,'xtick',[0:4:48])
set(gca,'XTickLabel',{'18','22','26','30','34','38','42','46','50','54','58','62','66'})
ylabel('Estimated CT')
set(gca,'ytick',[0:4:48])
set(gca,'YTickLabel',{'18','22','26','30','34','38','42','46','50','54','58','62','66'})



% following figure plots all the likelihood curves at all 12 pairs of times
figure()
hold on
for organ=1:8;
    for j =[2,4,6,8,10,12];
        hold on
        subplot(2,3,j/2)
        plot(linspace(0,46,769),reshape(M2(organ,:,j),769,1),'b');
    end
    for j =2:2:12;
        hold on
        subplot(2,3,j/2)
        plot(linspace(0,46,769),reshape(M2(organ,:,j+12),769,1),'b');
        %q=area(linspace(0,46,769),reshape(M2(organ,:,j+12),769,1));
        %alpha(q,0.2)
    end
end
for j =2:2:12;
    hold on
    subplot(2,3,j/2)
    plot([(j-1)*2,(j-1)*2],[0,7],'r-','Linewidth',3)
    plot([(j+11)*2,(j+11)*2],[0,7],'r-','Linewidth',3)
    axis([0 46 0 max(max(max(M2(:,:,[j j+12]))))])
    set(gca,'xtick',[0:4:48])
    xtickangle(90)
    set(gca,'XTickLabel',{'18','22','26','30','34','38','42','46','50','54','58','62','66'})
    xlabel('CT')
    ylabel('Likelihood')
end







