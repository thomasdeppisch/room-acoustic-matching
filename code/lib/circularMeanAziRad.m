function circMeanRad = circularMeanAziRad(aziRad)
% Calculate the circular mean, i.e., the mean of angles on a circle.
% https://en.wikipedia.org/wiki/Circular_mean

circMeanRad = atan2(sum(sin(aziRad)), sum(cos(aziRad)));