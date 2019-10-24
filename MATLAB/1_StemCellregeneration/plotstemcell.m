function plotstemcell()
A=load('output/md-ntpx.dat');
subplot(4,2,3);
loglog(A(:,1)/24,A(:,2));
xlim([1e0 1e4]);
end