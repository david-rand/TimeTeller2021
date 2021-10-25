function likeli_figs_against_Theta(Thetas, Liks, lowest, highest)

[sTheta, I] = sort(Thetas);
J = find(and(sTheta > lowest,sTheta < highest));
% I(J) are the indices of the individs with Theta < 0.05
figure
demoC(Liks(:,I(J)));
return

