function Yq = FastInterp(X, Y, Xq, inc)
%FASTINTERP is basically a faster version of interp1 with less features.
%Inputs are X, Y, X_query, and the increment in which the x data is
%provided in.

%Low and high indexes of query
indexFractional = (Xq - X(1)) / inc;
indexLow = floor(indexFractional) + 1;
indexHigh = ceil(indexFractional) + 1;

len_X = length(X);

% Make sure things don't go out of bounds
indexIsTooHigh = indexLow > len_X - 1;
indexLow(indexIsTooHigh) = len_X - 1;

indexIsTooLow = indexLow < 1;
indexLow(indexIsTooLow) = 1;

% Make sure things don't go out of bounds
indexIsTooHigh = indexHigh > len_X;
indexHigh(indexIsTooHigh) = len_X;

indexIsTooLow = indexLow < 2;
indexHigh(indexIsTooLow) = 2;

%Values at low and high indicies
XLow = X(indexLow);
XHigh = X(indexHigh);
YLow = Y(indexLow);
YHigh = Y(indexHigh);

%Interpolated value
Yq = (YHigh-YLow)/(XHigh-XLow)*(Xq-XLow)+YLow;

%Counters NaN when indicies are identical. 
if(indexLow==indexHigh)
    Yq = Y(indexLow);
end

end