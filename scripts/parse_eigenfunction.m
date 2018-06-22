function [Ns, nus, plot_x, plot_y] = parse_eigenfunction(path)
  sep = " ";
  Ns = [];
  nus = [];
  plot_x = [];
  plot_y = [];

  line = 0;
  header = dlmread(path, sep, [line, 0, line, 3]);
  while (!isempty(header))
    Ns = [Ns header(1)];
    nus = [nus header(2)];
    np = header(3);

    data = dlmread(path, sep, [line + 1, 0, line + np, 1]);
    plot_x = [plot_x plot(:, 1)];
    plot_y = [plot_y plot(:, 2)];

    line = line + np + 1;
    header = dlmread(path, sep, [line, 0, line, 3]);
  endwhile
endfunction
