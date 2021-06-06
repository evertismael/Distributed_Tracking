clear; clc;
tr = @(x,b,c)b/(1+exp(-x))-c;
b = 1;

syms x
expr = @(b,c) int(tr(x,b,c),x);

expr(1,0)
x=1;
subs(expr(1,0))
