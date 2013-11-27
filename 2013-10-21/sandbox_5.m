%%%%% Assignment Week 5 Sandbox %%%%%%


%% filter example
NT = 10;
b = [0.9 0.8];
a = [ 1 0.5];
x = zeros(NT,1);
x(1) = 1;
y = filter(b,a,x);
stem(y);