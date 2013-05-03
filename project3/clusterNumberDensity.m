% n(s,p) is the probability for any site to be a part of a cluster of size
% s
L = 400; 
N = 30;
n = 1:N;
pc = 0.59275;
pabove = pc + n/(180+N);
pbelow = pc - n(N:-1:1)/(180+N);
p = [pbelow, pc, pabove]
l = round(log2(L*L));

binStarts = unique(round(logspace(0,log10(L*L),100)));
l = length(binStarts)

results = zeros(l, length(p));
areas = zeros(L*L,1);
nExperiments = 100; 
%p = pc;
for k = 1:length(p)
    k
    areas = zeros(L*L,1);
    for i=1:nExperiments
        r = rand(L,L);

        z = r<p(k); % This generates the binary array
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
    binStarts = unique(round(logspace(0,log10(L*L),100)));
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
%figure
%[B,P] = meshgrid(bins,p);
%surf(B,P,results')
l = {};
for i=1:(N+1)
   loglog(binStarts,results(:,i))
   hold all;
   %l(i) = sprintf('p = %.2f', p(i));
   l{i} = [sprintf('p = %.4f', p(i))];
end
xlabel('s','Interpreter','LaTex')
ylabel('$n(s,p,L= 400)$','Interpreter','LaTex')
legend(l)

figure;
l = {};
for i=1:(N+1)
   loglog(binStarts,results(:,N+i))
   hold all;
   l{i} = [sprintf('p = %.4f', p(N+i))];
end
legend(l)
xlabel('s','Interpreter','LaTex')
ylabel('n(s,p,L=400)','Interpreter','LaTex')
%legend(l)

line = ones(length(binStarts),1)/2;
sxsibelow = zeros(1,N);
figure;
l = {};
tau = 1.89
for i=1:(N)
   rel = results(:,i)./binStarts(:).^(-tau);
   rel(isinf(rel)) = 0;
   loglog(binStarts,rel)
   hold all;
   %l(i) = sprintf('p = %.2f', p(i));
   [m,index] = max(rel);
   line = ones(length(binStarts),1)*max(rel)/2;
   loglog(binStarts,line)
   [xpos, ypos] = intersections(binStarts, rel, binStarts, line)
   ypos = ypos(xpos>binStarts(index));
   xpos = xpos(xpos>binStarts(index));
   
   i
   plot(xpos(length(xpos)), ypos(length(xpos)),'*')
   sxsibelow(i) = xpos(length(xpos));
   l{i} = [sprintf('p = %.4f', p(i))];
end
loglog(binStarts,line)
l{N+1} = [sprintf('s_xsi cutoff')];
xlabel('s','Interpreter','LaTex')
ylabel('$n(s,p,L= 400)/n(s,p_c, L = 400)$','Interpreter','LaTex')
legend(l)


sxsiabove = zeros(1,N);
figure;
l = {};
for i=1:(N)
   rel = results(:,N+1+i)./binStarts(:).^(-tau);
   rel(isinf(rel)) = 0;
   loglog(binStarts,rel)
   hold all;
   [m,index] = max(rel);
   line = ones(length(binStarts),1)*max(rel)/2;
%    ypos = ypos(xpos>binStarts(index));
%    xpos = xpos(xpos>binStarts(index));
    [xpos, ypos] = intersections(binStarts, rel, binStarts, line);
   %xpos = xpos(xpos > binStarts(index));
   
   i
   plot(xpos(length(xpos)), ypos(length(ypos)),'*')
   loglog(binStarts,line)
   %l(i) = sprintf('p = %.2f', p(i));
   l{i} = [sprintf('p = %.4f', p(i))];
   
   
   sxsiabove(i) = xpos(length(xpos));
end
loglog(binStarts, line)
l{N+1} = [sprintf('s_{xsi} cutoff')];
xlabel('s','Interpreter','LaTex')
ylabel('$n(s,p,L= 400)/n(s,p_c, L = 400)$','Interpreter','LaTex')
legend(l)

figure(10);
e1 = plot(pabove, sxsiabove,'o');
hold on
plot(pbelow, sxsibelow,'o')

figure;
logpabove = log10(pabove-pc);
logpbelow = log10(pc-pbelow);
logsabove = log10(sxsiabove);
logsbelow = log10(sxsibelow);
plot(logpabove, logsabove);
hold on
plot(logpbelow, logsbelow)

a = polyfit(logpabove, logsabove,1);
b = polyfit(logpbelow, logsbelow,1);
sigma = -1/mean([a(1), b(1)]);
propfactor = 10^mean([a(2), b(2)]);
p = linspace(0.0036,0.1,1000);
figure(10);
pa = pc + p;
pb = pc - p;
y = propfactor*p.^(-1/sigma);
e2 = plot(pa,y,'r');
hold on
plot(pb, (y),'r')
legend([e1,e2],'experimental data', sprintf('best fit curve, sigma = %.3f', sigma), 'Interpreter', 'LaTeX', 'FontSize', 14)

b = xlabel('p'); set(b, 'FontSize', 14)
a = ylabel('$s_{\xi}$','Interpreter','LaTeX');set(a,'FontSize',14)
