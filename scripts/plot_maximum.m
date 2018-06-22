function plot_maximum(Ns, plot_y)
  figure;
  scatter (Ns, max (plot_y, [], 1), 200, 'filled')
  axis ([Ns(1), Ns(end), 0, max(plot_y(:))])
endfunction
