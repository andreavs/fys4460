close all;
clear all;

%
% exwalk . m
%A.1. PERCOLATION 173
% Example of use of the walk routine
% Generate spanning cluster (l - r spanning )
tic;
a = toc;
lxs = [32,64,128,256,512,1024];
lys = [32,64,128,256,512,1024];
pc = 0.59275;
len = 0.3;
ps = linspace(0,len,10);
ps = [pc+ps];
M_SC = zeros(length(ps),length(lxs));
allps = zeros(length(ps),length(lxs));
ncount = 0;
perc = [];
nExperiments = 100;
for q=1:length(lxs)
    lx= lxs(q)
    ly = lys(q);
    ps = linspace(0,len,10);
    ps = [pc+ps];
    len = len/1.75;
    allps(:,q) = ps;
    for w=1:length(ps)
        p = ps(w)
        for j=1:nExperiments

            perc = [];
            ncount = 0;
            while ( isempty(perc))
                ncount = ncount + 1;
                if ( ncount >1000)
                    return
                end
                z = rand ( lx , ly ) < p ;
                [ lw , num ]= bwlabel (z ,4);
                perc_x = intersect ( lw (1 ,:) , lw ( lx ,:));
                perc = perc_x(perc_x ~= 0);
            end
            perc;
            % s = regionprops ( lw , 'Area' );
            % clusterareas = cat (1 , s.Area );
            % percarea = clusterareas(perc(1));
            % i = find ( clusterareas == percarea );
            zz = lw == perc(1) ;
            % zz now contains the spanning cluster
            %imagesc ( zz ); % Display spanning cluster
            % Run walk on this cluster
            [l , r ] = walk ( zz );
            zzz = l .* r ; % Find points where both l and r are non - zero
            zadd = zz + zzz ;
            M_SC(w,q) = M_SC(w,q) + sum(sum(zzz > 0));
    %         subplot (3 ,2 ,1) , imagesc ( zz );
    %         xlabel('a) The percolating cluster')
    %         subplot (3 ,2 ,2) , imagesc ( zadd );
    %         xlabel('b) The singly connected bonds highlighted')
    %         subplot (3 ,2 ,3) , imagesc ( zzz >0);
    %         xlabel('c) The singly connected bonds')
    %         subplot (3 ,2 ,4) , imagesc ( l +r >0);
    %         xlabel('d) The union of the paths of l/r walker')
    %         subplot (3 ,2 ,5) , imagesc (l);
    %         xlabel('e) the left walker')
    %         subplot (3 ,2 ,6) , imagesc (r);
    %         xlabel('f) the right walker')
        end
        M_SC(w,q) = M_SC(w,q)/(nExperiments*lx*lx);

    end
end
M_SC
% b = toc;
% b-a
% 
% l = log10(lxs);
% msc = log10(M_SC);
% c = polyfit(l,msc,1);
% y = figure(1);
% plot(l,msc)
% t = linspace(l(1),l(length(l)),1001);
% f = c(1)*t + c(2);
% hold on
% plot(t,f)
% xlabel('log10(L)')
% ylabel('log10(M_SC)')
% legend('M_SC experimental', sprintf('best fit, D_SC = %.4f',c(1)));
% 
% saveas(y,'test.fig')
l = {};
for i=1:length(lxs)
    figure(1);
    plot((allps(:,i)-pc)*lxs(i)^(3/4), lxs(i)^(1.25)*M_SC(:,i))
    %plot(log10(ps-pc),log10(M_SC(:,i)) + (3/4)*log10(lxs(i))  )
    hold all
    figure(2)
    plot((ps-pc), M_SC(:,i))
    hold all
    l{i} = [sprintf('L = %d', lxs(i))];
end
figure(1)
legend(l)
xlabel('(p-p_c)L^{D_{SC}}')
ylabel('L^{d-D_{SC}}P(p,L)')
figure(2)
xlabel('p-p_c')
ylabel('P(p,L)')
legend(l)