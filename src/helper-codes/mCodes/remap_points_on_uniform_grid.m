function f = remap_points_on_uniform_grid(x,y,z,f)

  % generate meshgrid
  [X, Y, Z] = meshgrid(x,y,z);
  % permute velocity data to match with meshgrid
  f_perm = permute(f,[2 1 3]);
  % calculate number of points required if grid has uniform spacing
  N = (max(y)-min(y))/mean(diff(y)) + 1;
  % generate y-co-ordinates with uniform grid spacing
  y_2 = linspace(min(y), max(y), N);
  % create meshgrid for querying data
  [Xq,Yq,Zq] = meshgrid(x,y_2,z);
  % get data at query points using interp3
  f_perm_q = interp3(X,Y,Z,f_perm,Xq,Yq,Zq);
  % undo permutation and copy to original array
  f = permute(f_perm_q,[2 1 3]);
  
end