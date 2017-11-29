function valid = checkKKT(alpha, T, f, idxl, idxr)

ekkt = 1e-6;
idx = find((idxr &   f.*T < 1-ekkt) | ...
           (idxl & f.*T > 1 + ekkt) | ...
           (alpha > ekkt & alpha < 1-ekkt & abs(f.*T-1) > ekkt) | ...
            (idxl & alpha < 1-ekkt) | ...
            (idxr & alpha > ekkt),1);

if (~isempty(idx) || abs(sum(T.*alpha)) > ekkt),
    disp(num2str(idx));
    [T(idx).*f(idx) alpha(idx)]
    valid = false;
else
    valid = true;
end