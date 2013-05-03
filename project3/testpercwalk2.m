%
% testpercwalk . m
%
% Generate spanning cluster (l - r spanning )
tic;
a = toc;
nsteplist = 2.^(6:15); %1:1000:1e4;
N = length(nsteplist);
nExperiments = 1000;

pc = 0.59275;
Np = 7;
ps = pc + linspace(0,0.1,Np);
Rsquared = zeros(Np,N);
for k = 1:Np
    p = ps(k)
    for i=1:N
        i
        nstep = nsteplist(i);
        for j=1:nExperiments
            lx = 400;
            ly = 400;
            nnstep = nstep + 1;
            ncount = 0;
            perc = [];
            while ( size ( perc ,1)==0)
                ncount = ncount + 1;
                if ( ncount >1000)
                    return
                end
                z = rand ( lx , ly ) < p ;
                [ lw , num ]= bwlabel (z ,4);
                perc_x = intersect ( lw (1 ,:) , lw ( lx ,:));
                perc = perc_x(perc_x ~= 0);
            end
            % s = regionprops ( lw , 'Area' );
            % clusterareas = cat (1 , s . Area );
            % maxarea = max ( clusterareas );
            % i = find ( clusterareas == maxarea );
            zz = (lw == perc(1));
            % zz now contains the spanning cluster
            %imagesc ( zz ) , axis equal , axis tight
            rz = 1.0* zz ;
            r = rand ( nnstep ,1);
            [w , n ] = percwalk ( rz ,r ,0);
            x = w (1 ,:);
            y = w (2 ,:);
            Rsquared(k,i) = Rsquared(k,i) + (x(length(x)) - x(1))^2 + (y(length(y))-y(1))^2;
        end
        Rsquared(k,i) = Rsquared(k,i)/nExperiments;
    end
end


b = toc; 
b-a
ns = log10(nsteplist);
rs = log10(Rsquared);
l = {};
D = zeros(1,Np);
for i=1:Np
    plot(ns,rs(i,:))
    l{2*i-1} = [sprintf('Rsquared experimental, p = %.03f', ps(i))];
    a = polyfit(ns,rs(i,:),1);
    t = linspace(ns(1), ns(N),100);
    f = a(1)*t + a(2);
    D(i) = 2/a(1);
    hold all;
    plot(t,f)
    l{2*i} = [sprintf('best fit, exponent = %.04f',a(1))];
end
xlabel('log10(N)')
% legend('Rsquared experimental', sprintf('best fit, exponent = %.04f',a(1)))
ylabel('log10(<R^2>)')
legend(l)
% hold on , plot (y , x );
% hold off

nu = 4/3; beta = 0.13;
mu = (D-2)*(nu -beta/2);
figure;
plot(log10(ps-pc), log10(mu))
figure()

rs = (Rsquared);
l = {};
for i=1:Np
    ns = log10(nsteplist*(pc-ps(i))^(-(beta-2*nu-mu(i))));
    plot(ns,rs(i,:))
    l{i} = [sprintf('Rsquared experimental, p = %.03f', ps(i))];
    hold all
end
xlabel('(N)')
ylabel('<R^2>')
legend(l)

