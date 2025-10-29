function [A, terms] = polynomialDesignMatrix(xy, maxDegree)
%POLYNOMIALDESIGNMATRIX Build a 2D polynomial design matrix.
%   [A, TERMS] = POLYNOMIALDESIGNMATRIX(XY, MAXDEGREE) returns the design
%   matrix for all bivariate monomials with total degree up to MAXDEGREE.
%   XY is an N-by-2 matrix with coordinates in [0, 1]^2. TERMS is a
%   MAXTERMS-by-2 matrix storing the exponent pairs [px, py] used in the
%   monomials.

arguments
    xy (:, 2) double
    maxDegree (1, 1) double {mustBeInteger, mustBeNonnegative}
end

nPts = size(xy, 1);
terms = generate_terms(maxDegree);
A = ones(nPts, size(terms, 1));

x = xy(:, 1);
y = xy(:, 2);

for i = 1:size(terms, 1)
    px = terms(i, 1);
    py = terms(i, 2);
    if px == 0 && py == 0
        A(:, i) = 1;
    else
        A(:, i) = (x .^ px) .* (y .^ py);
    end
end

end

function terms = generate_terms(maxDegree)
terms = [];
for px = 0:maxDegree
    for py = 0:(maxDegree - px)
        terms(end + 1, :) = [px, py]; %#ok<AGROW>
    end
end
end
