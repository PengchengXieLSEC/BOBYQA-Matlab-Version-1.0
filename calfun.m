function [F] = calfun (N, X)
  F = 0.0e0;
  for I = 4:2:N
    for J = 2:2:(I - 2)
      TEMP = (X(I - 1) - X(J - 1)) ^ 2 + (X(I) - X(J)) ^ 2;
      TEMP = max(TEMP, 1.0e-6);
      F = F + 1 / sqrt(TEMP);
    end
  end
end
