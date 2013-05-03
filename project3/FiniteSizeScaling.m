close all;
clear all;

Ls = [25,50,100,200,400,800];
N = length(Ls);
betas = zeros(1,10);
percolation = linspace(0,1,1001);
nExperiments = 10;
pi = zeros(1,N);

tic;
%p = percolation(900);
pc = 0.59275;

a = toc;
c = 0;
tic
b = toc
targets = [0.3,0.8];
Nt = length(targets)
results = zeros(Nt,N);
results2 = zeros(Nt,N);
l1 = {};
l2 = {};
for l = 1:2
    target = targets(l)
    
    %matlabpool(4);
    for m=1:N
        L = round(Ls(m));

        step = 0.25;
        p = 0.5;
        for k=1:15
            pi = 0;
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

                connectingClusters = union(ud,lr);
                %iscolumn(conectingClusters)
                connectingClusters = connectingClusters(connectingClusters ~= 0);
                pi = pi + length(connectingClusters);



            end
            pi = pi/nExperiments;
            if(pi < target)
               p = p + step;

            else
               p = p - step;
            end

            step = step/2;

        end
        L
        p
        pi
        results(l,m) = p;


    end
    figure(1)
    hold all
    plot(log10(Ls),results(l,:))
    title('x=0.8')
    xlabel('log10(L)')
    ylabel('p')
    pc = 0.59275;
    figure(2)
    Ls2 = log10(Ls);
    results2(l,:) = log10(abs(pc-results(l,:)));

    plot(Ls2, results2(l,:))
    hold all
    a = polyfit(Ls2, results2(l,:),1);
    t = linspace(Ls2(1), Ls2(N),100);
    f = a(1)*t + a(2);
    plot(t,f,'r')
    title('x=0.8')
    xlabel('log10(L)')
    ylabel('log10(p-pc)')
    l2{length(l2)+1} =  [sprintf('data set, x = %.03f', targets(l))];
    l2{length(l2)+1} = [sprintf('exponent = %.03f', a(1))];
    
end
legend(l2)
c = toc; 
c-b
figure(3);
nu = 4/3
l3 = {}
for i=1:N
    l3{i} = [sprintf('L = %d', Ls(i))];
    plot(Ls(i).^(1/nu)./(targets-pc), results(:,i))
    hold all
end
legend(l3)
