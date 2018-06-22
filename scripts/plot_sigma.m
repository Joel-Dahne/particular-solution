function plot_sigma(nus, sigmas)
  figure;
  hold on;
  for i=1:columns (nus)
    scatter(nus(:, i), sigmas(:, i), 'filled')
  endfor
endfunction
