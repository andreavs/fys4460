function [ beta ] = curveFitBeta( P, pc )
p = linspace(0,1,length(P))';

size(p)
size(P)
% p
model = @(beta) (p>pc).*(p-pc).^beta;
% beta = [0.18, 0.19]
% a = lsqcurvefit(@(beta) model(beta),beta,p,P)
% plot(p,model(0.2),'r')

nLargerThanPc = sum(p>pc) - 1;
P2 = P(length(P)-nLargerThanPc:length(P));
p2 = p(length(p)-nLargerThanPc:length(p)) - pc;
logp = log10(p2);
logP = log10(P2);

beta = polyfit(logp,logP,1);
factor = 10^beta(2);
beta = beta(1)

plot(p,P)
hold('on')
plot(p, factor*model(beta),'r')
figure()
% plot(logp, logP)
% figure()
% plot(p,P)
end

