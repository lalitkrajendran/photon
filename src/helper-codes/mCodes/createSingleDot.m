%% This code creates a single white dot on a black background

N = 1024;
A = zeros(N,N);

y = 768;
x = 256;

r = 4; % radius of the dot

A = MidpointCircle(A,r,x,y,255);

imshow(A)

fID = fopen('Dot_TopRight.bin','wb');
fwrite(fID,single(A'),'uint32');

