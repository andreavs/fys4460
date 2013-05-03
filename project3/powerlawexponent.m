nPoints = 1e6;
z = rand(nPoints,1).^(-3+1);

% crude atempt
% bar(hist(z) ./ sum(hist(z)))
% figure()

% better attempt
n = 51;
bins = zeros(n,1);
binStarts = logspace(0,10,n+1);
for i=1:length(bins)
    tester = (z>binStarts(i)).*(z<binStarts(i+1));
    bins(i) = sum(tester);
    binSize = (binStarts(i+1)-binStarts(i));
    bins(i) = bins(i)/binSize;
    
end  
bins = bins/nPoints;
starts = log10(binStarts(1:length(bins)));
log10bins = log10(bins)';
plot(starts,bins)
figure()
plot(starts,log10bins)
alpha = polyfit(starts,log10bins,1)


% culumative distribution
n = 1001;
bins = zeros(n,1);
starts = logspace(0,10,n)';
for i=1:n
    bins(i) = sum(z<starts(i));
    
end
bins = bins/nPoints;
figure()
plot(log10(starts),bins)
figure()
plot(log10(starts),log10(1-bins))
alphaplusone = polyfit(log10(starts),log10(1-bins),1)




density = zeros(n-1,1);
density(:) = bins(2:n) - bins(1:n-1);

binStarts = binStarts(1:n-1)';
density = density./binStarts;
figure()
plot(log10(binStarts),density);
alpha = polyfit(log10(binStarts),log10(density),1)
