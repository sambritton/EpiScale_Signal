function [length] = polyperi(X,Y)
edge_length = sqrt( ([X(2:end); X(1)]-X).^2 + ([Y(2:end); Y(1)]-Y).^2 );
length = sum(edge_length);