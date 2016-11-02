clear

r = 30e-3; % m
beta = 1.122*100; % 1 /m
k = linspace(0.1, 10, 1000)/1000; % m

plot(k*1000,exp(beta*(r-sqrt(r^2+k.^2))));
xlabel('Beam Separation (mm)')
ylabel('Energy from beam 2 / Energy from beam 1')

title('Energy of Interference from Second Beam At Beam 1')