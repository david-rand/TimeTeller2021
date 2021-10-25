function mouse_plot(Th2)

% This makes various plots from the output of the leave_1_tissue_out
% program. The first plot is of the likelihoods at six time points. The
% second plots the estimated time against the real time and the third plots
% a larger number of liklihood curves. Finally, its estimates the mean
% errors.

% History of modifications: Created from mouse_training_TimeTeller
% 26/12/19. Cleaned up 27/12/19 and used to produce figs for paper.

actuals = linspace(0,46,24)';

Likelis = Th2.Likelis_all;
num_finetimes = size(Th2.Likelis_all,2);

a=Th2.D_Ts;
a=a/769*48; %converts indices into times
%%
M2=Th2.Likelis_all;
%%
% following figure plots all the likelihood curves at all 12 pairs of times
figure()
hold on
% plot likelihood curves
for organ=1:8;
    for j =1:12;
        hold on
        subplot(4,3,j)
        %p=area(linspace(0,46,769),reshape(M2(organ,:,j),769,1));
        plot(linspace(0,46,num_finetimes),reshape(M2(organ,:,j),num_finetimes,1),'b');
    end
    for j =1:12;
        hold on
        subplot(4,3,j)
        plot(linspace(0,46,num_finetimes),reshape(M2(organ,:,j+12),num_finetimes,1),'b');
        %q=area(linspace(0,46,769),reshape(M2(organ,:,j+12),769,1));
        %alpha(q,0.2)
    end
end
%%
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
%%

% following figure plots estimated time T against actual time
figure()
% hold on
% for t = 1:24;
%     for organ= 1:8;
%         hold on
%         plot(actuals(t),a(t,organ),'b*','Markersize',6,'Linewidth',1)
%     end
% end
plot([-1 49], [-1 49],'k-')
hold on
plot([-1 25], [23 49],'k-')
plot([23 49], [-1 25],'k-')
plot([-1 1], [47 49],'k-')
plot( [47 49],[-1 1],'k-')
groups=repmat(1:8,24,1);
ractuals=repmat(actuals,1,8);
cm=jet(8);
gscatter(ractuals(:),a(:),groups(:),cm,'*',12,'on')
grid on
box on
axis([ -1 48 -1 48])
xlabel('Real CT')
set(gca,'xtick',[0:4:48])
set(gca,'XTickLabel',{'18','22','26','30','34','38','42','46','50','54','58','62','66'})
ylabel('Estimated CT')
set(gca,'ytick',[0:4:48])
set(gca,'YTickLabel',{'18','22','26','30','34','38','42','46','50','54','58','62','66'})
set(gcf,'Units','inches');
%
% uncoment this to output pdf figure - without the PaperPosition and
% PaperSize stuff Matlab has problems with this figure and pdf
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters epsFig
% The first two lines measure the size of your figure (in inches). The next line configures 
% the print paper size to fit the figure size. The last line uses the print command 
% and exports a vector pdf document as the output.

%%

% following figure plots all the likelihood curves at all 12 pairs of times
figure()
hold on
for organ=1:8;
    for j =[2,4,6,8,10,12];
        hold on
        subplot(2,3,j/2)
        plot(linspace(0,46,769),reshape(M2(organ,:,j),num_finetimes,1),'b');
    end
    for j =2:2:12;
        hold on
        subplot(2,3,j/2)
        plot(linspace(0,46,769),reshape(M2(organ,:,j+12),num_finetimes,1),'b');
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
%%
% Now we estimate the errors
k=[-48 -24 0 24 48];
for t = 1:24;
    for organ= 1:8;
        err(t,organ)=actuals(t)-a(t,organ);
    end; 
end
% The following deals with the fact that a may differ by actuals by
% soemthing close to a multiple of 24
for t = 1:24;
    for organ=1:8;
        [m,I]=min(abs(k-err(t,organ)));
        err2(t,organ)= k(I)-err(t,organ);
    end;
end
disp('mean abs deviation from actual in tissues')
mean(abs(err2))
disp('mean abs error')
mean(mean(abs(err2)))
disp('mean abs error from mean')
mean(abs(err2-mean(err2)))
disp('mean Thetas in tissues')
mean(Th2.D_Thetas)
%%
% Now a boxplot of the theta distribution fpor each tissue.
groups=1:8;
groups=repmat(groups,24,1);
figure
boxplot(Th2.D_Thetas(:),groups(:));
return