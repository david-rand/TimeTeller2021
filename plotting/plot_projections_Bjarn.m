function plot_projections_Bjarn(excluded_probes)

% This programme plots the data in 3d using the principal components from
% time 1 to 6.
% The input is a list of excluded probes. It is not necessary to exclude
% any here but when analysing the REMAGUS data we remove 15 as this has
% anomolous behaviour. See SI.

load('../data/frmaOM')
load('../data/Inx16_Clean')
%load('Inx16')
load('../data/Probes_string')
% datastr=load(data_str)
% names = fieldnames(datastr);

Inx = Inx16_Clean;
%Inx = Inx16;
if ~isempty(excluded_probes)
    Inx(excluded_probes)=[];
    % Inx([15])=[];
end
times = {'am8','midday','pm4','pm8','midnight','am4'};
data_times = [0 4 8 12 16 20];

cm=hsv(6);
cm = [1 0 0
    0 1 0
    1 0.5 0
    0 0 1
    1 0.4 1
    .2 1 1];
%cm=jet(6);
%cm=hsv(6);

clockOMstrN = make_clock_OM_str(frmaOM,Inx,Probes_string,[1:10],times,'Affymetrix U133 2.0');
%pm=probabilty_model_construct(clockOMstrN,Inx,[1:10],times,'Affymetrix U133 2.0');
%probabilty_model_construct(clockOMstr,rhythmic_probeset,instances_list,timeset,data_times,pm_period,pm_finetime)
pm=probabilty_model_construct(clockOMstrN,Inx,[1:10],{'am8','midday','pm4','pm8','midnight','am4'},data_times,24,7.5);

%set(gcf,'renderer','Painters')
for individual = 1:1
	%figure()
    for i = 1:6; % 1:24 % this is time as well
        %subplot(2,3,i)
        figure
        set(gcf,'renderer','Painters')
        for tms =1:6;%[1 3 6 9 12 15 18 21 24];
            hold on
            z=pm.proj.(times{i})*pm.datastr.(times{tms});
            plot3(z(1,:),z(2,:),z(3,:),'o','Color',cm(tms,:),'MarkerFaceColor',cm(tms,:),'Markersize',9)
            %plot3(pm.Mu.(times{i})(1,:),pm.Mu.(times{i})(1,:)(2,:),pm.Mu.(times{i})(1,:)(3,:))
            %plot3(pm.Mu.(times{i})(1,:),pm.Mu.(times{i})(1,:)(2,:),pm.Mu.(times{i})(1,:)(3,:))
            if i==tms
                plot3(z(1,:),z(2,:),z(3,:),'k.')
            end
        end
        mu = pm.Mu.(times{i});
        plot3(mu(1,:),mu(2,:),mu(3,:),'LineWidth',2)
        %plot3(PROJECTION.test.(times{i})(1,:),PROJECTION.test.(times{i})(2,:),PROJECTION.test.(times{i})(3,:),'k+','Markersize',10,'Linewidth',2)        
        xlabel('lPC1')
        ylabel('lPC2')
        zlabel('lPC3')
        grid on
        title({'Projection with PCs of Time',i})
    end
    legend([times,'Original projection', 'Test individual' ])
end