clear
h = [-1 2 -2 1];
xcorrH = -xcorr(h,h)
convH  = conv(h,h)

a = [-1 2 -1];
b = [1 -2 1];

xCorrA = -conv(a,a)
xCorrB = -conv(b,b)

