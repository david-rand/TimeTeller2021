%%
Th_e = human_TT(esetfrmaBC,[15],-15);
%%
[ r,tms ] = insert_real_times( remaguswithclock, Th_e.Tclock, Th_e.D_Thetas );
%%
I=find(tms(:,3)<0.7);
tms2=tms(I,:);
JJ=find(tms2(:,2)>10);
tms3=tms2(JJ,:);
%%
len=length(tms2(:,2));
clear prs
th = 0.02:0.001:0.16;
for i = 1:length(th)
    J(i) = length(find(tms3(:,3)<th(i))); % number of people with Th > th(i) in correct timing group
    K(i) = length(find(tms2(:,3)<th(i))); % number of people with Th > th(i)
    prs(i)=J(i)/K(i);
end
plot(th,prs)
%%
if secondpeak
    Tl=T-0.25;
    Tr = T+0.25;
    if Tl < 0
        sec_criticals = criticals(find(and(criticals>=0,criticals=<Tr);
        sec_criticals = [sec_criticals criticals(find(and(criticals>=T+0.75,criticals=<1)];
    else if Tr >1
        sec_criticals = criticals(find(and(criticals>=0,criticals=<T-0.75);
        sec_criticals = [sec_criticals criticals(find(and(criticals>=T-0.25,criticals=<1)];
        else
            sec_criticals = criticals(find(and(criticals>=Tl,criticals=<Tr);
        end
    end
end
            