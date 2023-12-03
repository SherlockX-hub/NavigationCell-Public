%% function of the correlation between two maps;
% Calculates the correlation for a point in the autocorrelogram. It is
% using the Pearsons correlation method. The 2 maps are assumed to be of
% same size.
function Rxy = zeroLagCorrelation(map1,map2)

[numRows, numCols] = size(map1);

sumXY = 0;
sumX = 0;
sumY = 0;
sumX2 = 0;
sumY2 = 0;
NB = 0;
for r = 1:numRows
    for c = 1:numCols
        if ~isnan(map1(r,c)) && ~isnan(map2(r,c))
            NB = NB + 1;
            sumX = sumX + map1(r,c);
            sumY = sumY + map2(r,c);
            sumXY = sumXY + map1(r,c) * map2(r,c);
            sumX2 = sumX2 + map1(r,c)^2;
            sumY2 = sumY2 + map2(r,c)^2;
        end
    end
end

if NB >= 0
    sumx2 = sumX2 - sumX^2/NB;
    sumy2 = sumY2 - sumY^2/NB;
    sumxy = sumXY - sumX*sumY/NB;
    if (sumx2<=0 && sumy2>=0) || (sumx2>=0 && sumy2<=0)
        Rxy = NaN;
    else
        Rxy = sumxy/sqrt(sumx2*sumy2);
    end
else
    Rxy = NaN;
end
end