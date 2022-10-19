function xodd = roundodd(x)

xfloor = floor(x);
xceil = ceil(x);

if mod(xfloor, 2) == 0
    xodd = xceil;
else
    xodd = xfloor;
end