%
% exflow . m
%
clear all ; clf ;
tic;
ag = toc;
% First , find the backbone
% Generate spanning cluster (l - r spanning )
lxs = [32,64,128,256, 512, 1024];
lys = [32,64,128,256, 512, 1024];
N = length(lxs);
pc = 0.59275;
Np = 1;
ps = pc;% + linspace(0,0.25,Np);
results = zeros(Np,N);
nExperiments = 500;
DEres = zeros(Np,N);
SCres = zeros(Np,N);
BBres = zeros(Np,N);
conduct = zeros(Np,N);


for l=1:N
    lx = lxs(l)
    ly = lys(l);
    for j=1:Np
        p = ps(j)
        for k=1:nExperiments
            ncount = 0;
            perc = [];
            while ( size ( perc ,1)==0)
                ncount = ncount + 1;
                if ( ncount >100000)
                    return
                end
                z = rand ( lx , ly ) < p ;
                [ lw , num ]= bwlabel (z ,4);
                perc_x = intersect ( lw (1 ,:) , lw ( lx ,:));
                perc = perc_x(perc_x ~= 0);
            end
%             s = regionprops ( lw , 'Area' );
%             clusterareas = cat (1 , s.Area );
%             maxarea = max ( clusterareas );
%             i = find ( clusterareas == maxarea );
            zz = lw == perc(1);
            % zz now contains the spanning cluster
            % Transpose
            zzz = zz';
            % Generate bond lattice from this
            g = sitetobond ( zzz );
            % Generate conductivity matrix
            [ p1 c_eff ] = FIND_COND (g , lx , ly );
            conduct(j,l) = conduct(j,l) + c_eff;
            % Transform this onto a nx x ny lattice
            x = coltomat ( full ( p1 ) , lx , ly );
            P = x .* zzz ;
            g1 = g (: ,1);
            g2 = g (: ,2);
            z1 = coltomat ( g1 , lx , ly );
            z2 = coltomat ( g2 , lx , ly );
            % Plotting

            f2 = zeros ( lx , ly );
            for iy = 1: ly -1
                f2 (: , iy ) = ( P (: , iy ) - P (: , iy +1)).* z2 (: , iy );
            end
            f1 = zeros ( lx , ly );
            for ix = 1: lx -1
                f1 ( ix ,:) = ( P ( ix ,:) - P ( ix +1 ,:)).* z1 ( ix ,:);
            end
            % Find the sum of absolute fluxes into each site
            fn = zeros ( lx , ly );
            fn = fn + abs ( f1 );
            fn = fn + abs ( f2 );
            fn (: ,2: ly ) = fn (: ,2: ly ) + abs ( f2 (: ,1: ly -1));
            fn (: ,1) = fn (: ,1) + abs (( P (: ,1) - 1.0).*( zzz (: ,1)));
            fn (2: lx ,:) = fn (2: lx ,:) + abs ( f1 (1: lx -1 ,:));

            limit = 0.0001;
            zfn = fn > limit;

            %[l , r ] = walk ( zz );
            epsilon = 1e-5*max(max(fn));
            zzzz = (abs(fn-2*sum(fn(:,1))) < epsilon).*(fn>0);
            BBres(j,l) = BBres(j,l) + sum(sum(zfn));
            DEres(j,l) = DEres(j,l) + sum(sum(zzz-zfn));
            SCres(j,l) = SCres(j,l) + sum(sum(zzzz));
            
            
            zbb =  ( zzz + zfn + zzzz);
            zbb = zbb / max ( max ( zbb ));
        end
        BBres(j,l) = BBres(j,l)/(nExperiments);
        DEres(j,l) = DEres(j,l)/(nExperiments);
        SCres(j,l) = SCres(j,l)/(nExperiments);
        conduct(j,l) = conduct(j,l)/nExperiments;
    end
end

bb = log10(BBres);
dd = log10(DEres);
sc = log10(SCres);
cc = log10(conduct);
ll = log10(lxs);
l = {};
l2 = {};
l3 = {};
t = linspace(ll(1),ll(N),11);
for i=1:Np
    figure(1);
    plot(ll,bb(i,:));
    a = polyfit(ll,bb(i,:),1);
    f = a(1)*t + a(2);
    hold all
    plot(t,f)
    l{8*i-7} = [sprintf('M_{BB}, p = %.04f', ps(i))];
    l{8*i-6} = [sprintf('best fit, D = %.04f', a(1))];
    
    
    %figure(2)
    plot(ll,dd(i,:))
    a = polyfit(ll,dd(i,:),1);
    f = a(1)*t + a(2);
    hold all
    plot(t,f)
    l{8*i-5} = [sprintf('M_{DE}, p = %.04f', ps(i))];
    l{8*i-4} = [sprintf('best fit, D = %.04f', a(1))];
    
    %figure(3)
    plot(ll,sc(i,:))
    hold all
    a = polyfit(ll,sc(i,:),1);
    f = a(1)*t + a(2);
    plot(t,f)
    l{8*i-3} = [sprintf('M_{SC}, p = %.04f', ps(i))];
    l{8*i-2} = [sprintf('best fit, D = %.04f', a(1))];
    
    plot(ll,cc(i,:))
    hold all
    a = polyfit(ll,cc(i,:),1);
    f = a(1)*t + a(2);
    plot(t,f)
    l{8*i-1} = [sprintf('conductivitiy, p = %.04f', ps(i))];
    l{8*i} = [sprintf('best fit, zeta = %.04f', a(1))];
    i
    
end
figure(1)
legend(l)
% figure(2)
% legend(l2)
% figure(3)
% legend(l3)
% figure
% subplot (2 ,2 ,1) , imagesc ( zzz );
% title ( ' Spanning cluster ')
% axis equal
% axis tight
% subplot (2 ,2 ,2) , imagesc ( P );
% title ( ' Pressure ' );
% axis equal
% axis tight
% subplot (2 ,2 ,3) , imagesc ( fn );
% title ( 'Flux' );
% axis equal
% axis tight
% subplot (2 ,2 ,4) , imagesc ( zbb );
% title ( ' BB, SC and DE ' );
% axis equal
% axis tight
% 
% figure;
% zadd = zz + zzzz ;
% subplot (3 ,2 ,1) , imagesc ( zz );
% xlabel('a) The percolating cluster')
% subplot (3 ,2 ,2) , imagesc ( zadd );
% xlabel('b) The singly connected bonds highlighted')
% subplot (3 ,2 ,3) , imagesc ( zzzz >0);
% xlabel('c) The singly connected bonds')
% subplot (3 ,2 ,4) , imagesc ( l +r >0);
% xlabel('d) The union of the paths of l/r walker')
% subplot (3 ,2 ,5) , imagesc (l);
% xlabel('e) the left walker')
% subplot (3 ,2 ,6) , imagesc (r);
% xlabel('f) the right walker')

b = toc
b-ag