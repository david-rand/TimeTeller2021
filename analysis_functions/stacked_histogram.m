%%
d1 = Th.D_Thetas;                                           % Create Data
d2 = Th_pm.D_Thetas;                                           % Create Data
binrng = 0:0.05:1;                                                   % Create Bin Ranges
counts1 = histc(d1, binrng);                                    % Histogram For ?d1?
counts2 = histc(d2, binrng);                                    % Histogram For ?d2?
counts3 = counts1 + counts2;                                    % Histogram Sum ?d1?+?d2?
figure
bar(binrng, counts3, 'b')
hold on
bar(binrng, counts1, 'r')
hold off
legend('Training Data', 'REMAGUS')
%%