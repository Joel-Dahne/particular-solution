function [Ns, nus, sigmas] = parse_sigma(path)
  sep = " ";
  Ns = [];
  nus = [];
  sigmas = [];

  line = 0;
  header = dlmread(path, sep, [line, 0, line, 2]);
  while (!isempty(header))
    Ns = [Ns header(1)];
    np = header(2);

    data = dlmread(path, sep, [line + 1, 0, line + np, 1]);
    nus = [nus data(:, 1)];
    sigmas = [sigmas data(:, 2)];

    line = line + np + 1;
    header = dlmread(path, sep, [line, 0, line, 2]);
  endwhile
endfunction
