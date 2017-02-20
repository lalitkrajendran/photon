%% This code creates a dot card of white markers on a black background

N = 1024;
A = zeros(N,N);

del_x = 1024; % interval between dots
del_y = 1024; % interval between dots
r = 1; % radius of the dot

for i = del_x/2:del_x:N-del_x/2
    for j = del_y/2:del_y:N-del_y/2
        A = MidpointCircle(A,r,i,j,255);
%         A(i,j) = 255;
    end
end
imshow(A')