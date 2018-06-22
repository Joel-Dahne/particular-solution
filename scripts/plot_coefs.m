function plot_coefs(coefs)
  figure;
  hold on;
  for i=1:columns (coefs)
    scatter(1:length(coefs{i}), log(coefs{i}), 'filled')
  endfor
endfunction
