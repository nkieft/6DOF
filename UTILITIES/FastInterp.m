function Yq = FastInterp(X, Y, Xq, inc)
%FASTINTERP is basically a faster version of interp1 with less features.
%Inputs are X, Y, X_query, and the increment in which the x data is
%provided in.

%Low and high indexes of query
indexLow = floor(Xq/inc)+1;
indexHigh = ceil(Xq/inc)+1;

%Make sure things don't go out of bounds
if(indexLow > length(X))
    indexLow = length(X)-1;
end

%Make sure things don't go out of bounds
if(indexHigh > length(X))
    indexHigh = length(X);
end

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