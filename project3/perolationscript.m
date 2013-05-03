close all;
clear all;

Ls = logspace(2,4,10);
betas = zeros(1,10);
percolation = linspace(0,1,1001);
nExperiments = 50;
P = zeros(1001,1);

tic;
%p = percolation(900);
pc = 0.59275;

a = toc;
c = 0;
%matlabpool(4);
for m=1:10
    L = round(Ls(m))
    P = zeros(1001,1);
    for k=1:1001
        p = percolation(k);
        for i=1:nExperiments
            r = rand(L,L);

            z = r<p; % This generates the binary array
            [lw,num] = bwlabel(z,4);


            %img = label2rgb(lw,'jet','k','shuffle');
            %image(img);

            %s = regionprops(lw,'BoundingBox');
            %bbox = cat(1,s.BoundingBox);
            up = lw(1,:);
            down = lw(L,:);
            left = lw(:,1);
            right = lw(:,L);
            ud = intersect(up,down);
            lr = intersect(left,right);

            conectingClusters = union(ud,lr);
            %iscolumn(conectingClusters)
            l = length(conectingClusters);

            for j=1:l
                if(conectingClusters(j) ~= 0)
                    P(k) = P(k) + sum(sum(lw == conectingClusters(j)))/L^2;
                end
            end
            

        end

    end
    c = c/nExperiments;
    P = P/nExperiments;
    % plot(percolation,P)
    % figure()
    beta = curveFitBeta(P,pc)
    betas(m) = beta

end
c = c/nExperiments;
P = P/nExperiments;
b = toc;
plot(percolation,P);

b-a
plot(L,betas)