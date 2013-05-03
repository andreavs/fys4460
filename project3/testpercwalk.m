%
% testpercwalk . m
%
% Generate spanning cluster (l - r spanning )
tic;
a = toc;
nsteplist = 1e4;
N = length(nsteplist);
nExperiments = 1;
Rsquared = zeros(1,N);

    i = 1
    nstep = nsteplist(i);
    for j=1:nExperiments
        lx = 100;
        ly = 100;
        p = 0.59275;
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
        imagesc ( zz ) , axis equal , axis tight
        rz = 1.0* zz ;


        r = rand ( nnstep ,1);
        [w , n ] = percwalk ( rz ,r ,0);

        x = w (1 ,:);
        y = w (2 ,:);
        Rsquared(i) = Rsquared(i) + (x(length(x)) - x(1))^2 + (y(length(y))-y(1))^2;
    end
    Rsquared(i) = Rsquared(i)/nExperiments;

b = toc; 
b-a
ns = log10(nsteplist);
rs = log10(Rsquared);
% plot(ns,rs)
% a = polyfit(ns,rs,1);
% t = linspace(ns(1), ns(N),100);
% f = a(1)*t + a(2);
%hold all;
%plot(t,f)
%xlabel('log10(N)')
%legend('Rsquared experimental', sprintf('best fit, exponent = %.04f',a(1)))
%ylabel('<R^2>')
hold on , plot (y , x );
hold off