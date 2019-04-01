function [skewness] = tri_skewness(X,Y)
perimeter = polyperi(X,Y);
area      = polyarea(X,Y);
skewness  = abs(sqrt(3)/4*(perimeter/3).^2-area)./(sqrt(3)/4*(perimeter/3).^2);