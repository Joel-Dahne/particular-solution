function [Ns, nus, coefs] = parse_coefs(path)
  sep = " ";
  Ns = [];
  nus = [];
  coefs = {};

  line = 0;
  header = dlmread(path, sep, [line, 0, line, 2]);
  while (!isempty(header))
    Ns = [Ns header(1)];
    nus = [nus, header(2)];

    data = dlmread(path, sep, [line + 1, 0, line + header(1), 1]);
    coefs = {coefs{:} data};

    line = line + header(1) + 1;
    header = dlmread(path, sep, [line, 0, line, 2]);
  endwhile
endfunction
