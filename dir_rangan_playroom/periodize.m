function output = periodize(input,a,b);
output = a+mod(input-a,(b-a));
