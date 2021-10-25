function Bjarnason_timing_plotter_v2(Th,excluded_probes)
% plots the timing scatter plot of actual vs estimated

load('Probes_string')
load('Inx16_Clean')
Inx = Inx16_Clean;
%Inx = Inx16;
if ~isempty(excluded_probes)
    Inx(excluded_probes)=[];
    % Inx([15])=[];
end

probenames = Probes_string(Inx);
%just making the names easier to read in figures
probenames  = strrep(probenames , 'g', '');
probenames  =strrep(probenames , '_', ' ');
probenames  = strrep(probenames , 'at', '');
probenames  = strrep(probenames , 's', '');
probenames  = strrep(probenames , 'x', '');


num_finetimes = size(Th.Likelis_all,2);

% We use 16 probes here and the corresponding data from Bjarnason et al 
% is contained in bjarn16
load('Bjarn16');
clear clock_OM
clock_OM = bjarn16;ww


bmal=clock_OM(7,:); %data for  Bmal1

% use cosinor to get the phases of the probes for each individual
clear q a p allphases maxf minf antiphase
for q=1:length(Inx); % runs over probes
    for a = 1:10; % runs over individuals
        [p(q,a), allphases(q,a),maxf(q,a),antiphase(q,a), minf(q,a)]= cosinor(linspace(0,1,7),clock_OM(q,[(a-1)*6+1:6*a, (a-1)*6+1]),2*pi,0.05);
        if allphases(q,a)>24
            allphases(q,a)=allphases(q,a)-24;
        end
    end
end

a2 = Th.D_Ts/193 *24;
actuals = repmat(0:4:20,10,1)';
%actuals = [8 12 16 20 0 4]';

errors = actuals-a2; %difference estimated time and actual time
% must allow for fact that 24h error is not error
for t = 1:60;
    if errors(t)>12
        errors(t)=errors(t)-24;
    end
    if errors(t)<-12
        errors(t)=errors(t)+24;
    end
end

% merrors = mean(errors);

allphases(11,1)= allphases(11,1)-24;
allphases(12,1)= allphases(12,1)-24;
allphases(13,1)= allphases(13,1)-24;
% 
allphases(11,3)= allphases(11,3)-24;
allphases(12,3)= allphases(12,3)-24;
allphases(11,9)= allphases(11,9)-24;
allphases(14,1)= allphases(14,1)-24;
%  

% what folows plots the array of figures that show for each gene and
% individual the residual (x-axis) against the Cosinor phase of the gene
% (y-axis)
allphases(3,1)= allphases(3,1)-24;
cm=hsv(10);
figure()
hold on
for k=1:length(Inx); % running through probes
    subplot(4,4,k)
    plot(mean(errors),allphases(k,:),'*','Markersize',8,'Linewidth',1.5)
    lsline
    for indiv=1:10;
        hold on
        plot(mean(errors(:,indiv)),allphases(k,indiv),'*','Color',cm(indiv,:),'Markersize',8,'Linewidth',1.5)
    end

    [f gof]=fit(mean(errors)',allphases(k,1:10)','poly1');

    str=sprintf('f(x)= %.1f x+ %.1f R2= %.2f', f.p1 ,f.p2,gof.rsquare);
    hold on
    title({probenames{k}; str})
    plot([0 0],[min(allphases(k,:))-1 max(allphases(k,:))+1],'k-')
    plot([-4 3],[f.p2 f.p2],'k-')
    axis([-4 3 min(allphases(k,:))-1 max(allphases(k,:))+1])

    xlabel('Residual')
    ylabel('Cosinor phase')
    box on
    grid on
end
% end of plot

% Now we plot the likelihood curves for each timepoint in an array of
% figures
a=a2/193*24;
figure()
hold on
for individual=1:10;
for j =1:6;
    hold on
subplot(2,3,j)
plot(linspace(0,24,193),reshape(Mall(individual,:,j),193,1),'Color',cm(individual,:),'Linewidth',2);
end
end
legend('M15','M11','M12','M11','M9','F18','F13','F14','F5','F6')
for j =1:6;
    hold on
subplot(2,3,j)
plot([(j-1)*4,(j-1)*4],[0,6],'k-','Linewidth',2)
axis([0 24 0 max(max(Mall(:,:,j)))])
set(gca,'XTick' ,0:4:24)
xtickangle(90)
set(gca,'XTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
xlabel('Time')
ylabel('Likelihood')
end
subplot(2,3,1)
plot([(7-1)*4,(7-1)*4],[0,6],'k-','Linewidth',2)


% Now we plot actual time vs estimated time in a scatter plot with individuals color-coded 
actuals = linspace(0,20,6)';
figure()
hold on
for i=1:10;
    plot(actuals, a2(:,i),'*','Color',cm(i,:),'Markersize',8,'Linewidth',1.5)
end
plot([-1 25], [-1 25],'k-')
plot([-1 1], [23 25],'k-')
plot([23 25], [-1 1],'k-')
grid on
box on
axis([ -1 25 -1 25])
xlabel('Real Time')
set(gca,'XTick' ,0:4:24)
set(gca,'XTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
ylabel('Estimated Time')
set(gca,'YTick' ,0:4:24)
set(gca,'YTickLabel',{'8am','12day','4pm','8pm','12night','4am','8am'})
legend( 'M15','M11','M12','M11','M9','F18','F13','F14','F5','F6')
return