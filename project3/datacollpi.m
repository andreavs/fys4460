close all;
clear all;

Ls = [25,50,100,200,400,800,1600];
ps = [0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75];

N = length(Ls);

Np = 50;
ps = zeros(N,Np);
spread = 0.25;
pc = 0.59275;
for i=1:N
    p = linspace(pc-spread, pc+spread,Np);
    ps(i,:) = p;
    spread = spread/2;
end

nExperiments = 10000;


tic;
%p = percolation(900);

b = toc
results = zeros(Np,N);
results2 = zeros(2,N);
for l = 1:Np
    
    
    %matlabpool(4);
    for m=1:N
        p = ps(m,l)
        L = round(Ls(m));
        pi = 0;
        for k=1:nExperiments
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

            connectingClusters = union(ud,lr);
            %iscolumn(conectingClusters)
            connectingClusters = connectingClusters(connectingClusters ~= 0);
            pi = pi + min(1,length(connectingClusters));
            
       end
       pi = pi/nExperiments
       results(l,m) = pi;
       

    end

    
end
c = toc; 
c-b

nu = 4/3;
l = {};
for i=1:N
    plot((ps(i,:)-pc).*Ls(i).^(1/nu),results(:,i))
    hold all
    l{i} = [sprintf('L = %d', Ls(i))];
end
legend(l)
