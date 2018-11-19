x = [(1:5)/5 zeros(1,128-5)];
x = x(randperm(128));
y = x + 0.05*randn(1,128);
lambda = 0.2;
func = @(z) 0.5*sum((z-y).^2) + 0.5*lambda*sum((z.^2));
fplot(@(z) func(z))
[z,MIN] = fminbnd(func,-5,5);