function cdfplot(x)

sx=sort(x);
y=1:length(sx);
y=y/length(sx);
sx=[0 sx];
y=[0 y];
plot(y,sx);
xlim([0 1])
return