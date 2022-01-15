function out = interpolate(data, index)
  c = ceil(index);
  f = floor(index);
  fraction = c - f;
  out = (1 - fraction) * data(f,  :) + fraction * data(c, :);
end

