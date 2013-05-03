f = @(x) [2*(x(1)+1)*(x(2)+1)^2; 2*(x(2)+1)*(x(1)+1)^2];

first = [3;2];
g1 = f(first)
h = 0.001; 
hess1 = zeros(2);
g1plus = f([first(1)+h; first(2)])
g1minus = f([first(1)-h; first(2)])

hess1(1,1) = (g1plus(1) - g1minus(1))/(2*h);

g1plus = f([first(1); first(2)+h])
g1minus = f([first(1); first(2)-h])
hess1(2,2) = (g1plus(2) - g1minus(2))/(2*h);

hess1(1,2) = (g1plus(1) - g1minus(1))/(2*h);
g1plus = f([first(1)+h; first(2)])
g1minus = f([first(1)-h; first(2)])
hess1(2,1) = (g1plus(2) - g1minus(2))/(2*h);

B = hess1;

hess1
N = 20;


g1old = g1;
step = 1;
for i=1:N
    i
   
   p = -B^(-1)*g1
   first = first + step*p
   s = step*p
   g1 = f(first)
   y = g1-g1old
   1/(s'*s)*(y-B*s)*(s')
   B = B + 1/(s'*s)*(y-B*s)*s'
   g1old = g1;
    
    
end
