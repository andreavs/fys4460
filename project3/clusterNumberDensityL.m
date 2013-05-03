% n(s,p) is the probability for any site to be a part of a cluster of size
% s
Len = 2.^(4:11); 
N = 10;
n = 1:N;
pc = 0.59275;
pabove = pc + n/(80+N);
pbelow = pc - n(N:-1:1)/(80+N);
p = pc;
%l = round(log2(L*L));

binStarts = unique(round(logspace(0,log10(Len(length(Len))*Len(length(Len))),100)));
l = length(binStarts)

results = zeros(l, length(Len));
%areas = zeros(L*L,1);
nExperiments = 1000; 
%p = pc;
massresults = zeros(length(Len), nExperiments);
for k = 1:length(Len)
    k
    L = Len(k)
    areas = zeros(Len(length(Len))*Len(length(Len)),1);
    for i=1:nExperiments
        r = rand(L,L);

        z = r<p; % This generates the binary array
        [lw,num] = bwlabel(z,4);

        %img = label2rgb(lw,'jet','k','shuffle');
        %image(img);

        s = regionprops(lw,'Area');
        area = cat(1,s.Area);
        up = lw(1,:);
        down = lw(L,:);
        left = lw(:,1);
        right = lw(:,L);
        ud = intersect(up,down);
        lr = intersect(left,right);

        connectingClusters = sort(union(ud,lr));

        connectingClusters = connectingClusters(connectingClusters ~= 0);
        if(length(connectingClusters) > 0)
            massresults(k,i) = max(area(connectingClusters));
        end
        area(connectingClusters) = [];
        
%         counter = 0; % when we delete elements we must shift indices.
%         for j=1:length(connectingClusters)
%            if connectingClusters(j) ~= 0
%                index = connectingClusters(j) + counter;
%                area(index) = [];
%                counter = counter - 1; 
% 
%            end
%         end
        %areas
        c = hist(area,max(area)-min(area)+1);
        areas(min(area):max(area)) = areas(min(area):max(area)) + c';

        %area
        %areas
    end
    areas = areas/(nExperiments*L*L);%.*(1:L*L)';

    %n = 20
    %bins = zeros(n,1);

    %bins = zeros(1,1);
    %binStarts = zeros(1,1);
    binStarts = unique(round(logspace(0,log10(Len(length(Len))*Len(length(Len))),100)));
    l = length(binStarts);
    bins = zeros(1,l);
    for i=1:(l-1)
        binSize = binStarts(i+1) - binStarts(i);
        %areas(binStarts(i):(binStarts(i+1)-1));

        bins(i) = sum(areas(binStarts(i):(binStarts(i+1)-1)) )/binSize;
        
    end
    

%     binSize = 1;
%     counter = 1;
%     while binSize < (length(areas)/2)
%         %counter
%         bins(counter) = sum(areas(binSize:2*binSize-1))/binSize;
%         binStarts(counter) = binSize;
%         binSize = 2*binSize;
%         counter = counter+1;
% 
%     end
    results(:,k) = bins;
end
l = {}
for i=1:length(Len)
   loglog(binStarts,results(:,i))
   hold all;
   %l(i) = sprintf('p = %.2f', p(i));
   l{i} = [sprintf('L = %d', Len(i))];
end
xlabel('s')
ylabel('n(s,p=p_c,L)')
legend(l)



n = results(:,length(Len));
n = log10(n(1:70));
binStarts = log10(binStarts(1:70));
a = polyfit(binStarts,n',1)

massmeans = zeros(N,1);
for i=1:N
    b = massresults(i,:);
    b(b==0) = [];
    massmeans(i) = mean(b);
    
end
figure()
loglog(Len, massmeans)
loglog(Len, massmeans)
title('Mass of the percolating clusters')
xlabel('System length')
ylabel('Mass')
logmean = log10(massmeans);
logLens = log10(Len);

a = polyfit(logLens,logmean',1)

