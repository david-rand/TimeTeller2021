function xsrt=make_expression_str(probelist)

clock_OM = frmaOM(Inx16,:);
all_OM = frmaOM;
times = {'am8','midday','pm4','pm8','midnight','am4'};

for probe = probelist
    for tms = 1:6;
          for indivs=1:10;
            expressOMstr.(times{tms})(indivs) = all_OM(probe,6*(indivs-1)+tms);
          end
    end
end
return


function ls = test4rhythm(ilist)

for tms = 1:6;
      for indivs=1:10;
        expressOMstr.(times{tms})(indivs) = frmaOM(probe,6*(indivs-1)+tms);
      end
    end
end

count = 0;
for st=list
    for indiv=1:10
        for time = 1:6
            gene=expressOMstr.(times{time})(:,indiv);
            if 0
                s =std(genevect);
                me = mean(genevect);
                gene=(genevect(gnum)-me)/(s);
            else
                gene=genevect(gnum);
            end
            T = bjarnason.T(time,indiv);
            ind(count)=indiv;
            ls(1,count) = T+6;
            ls(1,count) = mod(ls(1,count),24);
            ls(2,count) = gene;
            count=count+1;
        end
    end
end



function index = get_indices
gene_names = {}
index=[];
strcmp(Probes_string{i}(1:3),'HLF');

for j=1:length(Probes_string)
    for i=1:length(Probes_string)
        if strcmp(Probes_string{i}(1:length(gene_names{j}), gene_names{j}(1:length(gene_names{j}))
            index = [index i];
        end
    end
end

% associate time T and Theta to each sample in regagous

% plot gene against T

% get curve 

% test curve for variance.
function [lst, nlst] = getlists(gnum,remagous,gtlt,theta)
lst=zeros(2,226);
nlst=zeros(2,226);
indivs=[]
if gtlt
    for indiv=1:226
        if (remegus.Theta(indiv))<theta)
            indivs=[indivs indiv];
        end
    end
else
     for indiv=1:226
        if (remegus.Theta(indiv))>theta)
            indivs=[indivs indiv];
        end
     end
end
    
for indiv=individs
    lst=zeros(2,226);
    T = remegus.Theta;
    express=remegus.expression(gnum,indiv);
    nexpress=remegus.nexpression(gnum,indiv);
    lst(1,indiv) = T+6;
    lst(2,indiv) = express;
    nlst(1,indiv) = T+6;
    nlst(2,indiv) = nexpress;
end
return

% ================
    function s=stdfitex(lst)
        
    f=fit(ls(1,:)',ls(2,:)','poly3');
    mi = min(ls(1,:));mx = max(ls(1,:));
    domain = mi:(mx-mi)/100:mx;
    s1=std(f(domain);
    s2=std(ls(2,:)-f(ls(1,:)))
    ratio = s1/s2
    return
    
    
    if remegus.Theta(indiv)>0.3
        genevect=remegus.expression(:,indiv);
        if 1
            s =std(genevect);
            me = mean(genevect);
            gene=(genevect(gnum)-me)/(s);
        else
            gene=genevect(gnum);
        end

        T = remegus.T(indiv)+6;
        lst2 = [lst2,[T, gene]'];
    end
end
scatter(lst2(1,:),lst2(2,:),'Color',[0 0.5 0.1])
hold off