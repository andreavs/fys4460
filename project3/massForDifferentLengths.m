N = 8;
massmeans = zeros(N,1);
for i=1:N
    b = massresults(i,:);
    b(b==0) = [];
    massmeans(i) = mean(b);
    
end

loglog(Len, massmeans)
title('Mass of the percolating clusters')
xlabel('System length')
ylabel('Mass')

logmean = log10(massmeans);
logLens = log10(Len);

a = polyfit(logLens,logmean',1)