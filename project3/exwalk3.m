% tic;
% a = toc;
%
% exwalk . m
%A.1. PERCOLATION 173
% Example of use of the walk routine
% Generate spanning cluster (l - r spanning )
lx =64;
ly = 64;
p = 0.59275;
p = 0.7
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
%s = regionprops ( lw , ' Area ' );
%clusterareas = cat (1 , s . Area );
%maxarea = max ( clusterareas );
%i = find ( clusterareas == maxarea );
zz = lw == perc(1) ;
% zz now contains the spanning cluster
imagesc ( zz ); % Display spanning cluster
% Run walk on this cluster
[l , r ] = walk ( zz );
zzz = l .* r ; % Find points where both l and r are non - zero
zadd = zz + zzz ;
subplot (2 ,3 ,1) , imagesc (zz );
subplot (2 ,3 ,2) , imagesc (zadd );
subplot (2 ,3 ,3) , imagesc (zzz >0);
subplot (2 ,3 ,4) , imagesc (l +r >0);
subplot (2 ,3 ,5) , imagesc (l);
subplot(2,3,6) , imagesc(r);

% b = toc;
% b-a