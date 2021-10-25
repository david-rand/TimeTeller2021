function [ r,tms ] = insert_real_times( rin, Ts, Thetas, Ztimes )
% This produces a table with the real times inserted and also outputs an
% array with the real times in the first column, the estimate times in
% the second and the Thetas in the third.
% Inputs
% rin a table that the real times will be inserted into
% Ts the estimated times to be inserted (column)
% Theta values state will also be inserted (column same size as Ts)
% Ztimes: the real tiems to be inserted

m = 0;
for i = 1:length(rin.GSM) % any colimn to get no of patients
	m = max(m,str2num(convertStringsToChars(rin.VarName1(i))));
end
% m gives maximum patient number

rtimes = -ones(m,1);
M=csvread('realtime.csv',1,0);
M(:,1)=int16(M(:,1)); % col 1 is the patient id

for i=1:size(M,1) rtimes(M(i,1))=M(i,2);end
% rtimes(patient number) is that patients sampling time
J=find(rtimes>0);
tms = zeros(length(J),5);

count = 0;
 for i = 1:length(rin.VarName1) % rin.VarName1 is patient id
     indx = str2num(convertStringsToChars(rin.VarName1(i)));
     t=rtimes(indx); % real time of sample for that patient
     if t>-1 % -1 if no real time available
        r.realtimes(i) = t;
        count = count + 1;
        tms(count,1)=t;
        tms(count,2)=Ts(i);
        tms(count,3)=Thetas(i);
        tms(count,4)=i;
        tms(count,5)=Ztimes(i);
     else
         r.realtimes(i) = NaN;
     end
     r.index(i)=i;
 end
end

