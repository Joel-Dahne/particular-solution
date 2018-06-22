function plot_eigenfunction(plot_x, plot_y)
  figure;
  hold on;
  for i=1:columns (plot_x)
    scatter(plot_x(:, i), plot_y(:, i), 'filled')
  endfor
endfunction
