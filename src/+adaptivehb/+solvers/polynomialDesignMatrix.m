function [A, terms] = polynomialDesignMatrix(xy, maxDegree)
%POLYNOMIALDESIGNMATRIX Build a 2D polynomial design matrix.
%   [A, TERMS] = POLYNOMIALDESIGNMATRIX(XY, MAXDEGREE) returns the design
%   matrix A for all bivariate monomials x^px * y^py whose total degree
%   px + py does not exceed MAXDEGREE.
%
%   Inputs:
%       XY        - N-by-2 matrix of coordinates (typically in [0,1]^2).
%       MAXDEGREE - non-negative integer setting the maximum total degree.
%
%   Outputs:
%       A     - N-by-M design matrix where M = (maxDegree+1)*(maxDegree+2)/2.
%               Each column corresponds to one monomial term.
%       TERMS - M-by-2 matrix of exponent pairs [px, py].
%
%   The terms are enumerated in graded lexicographic order.  For example,
%   with maxDegree = 2 the columns of A correspond to:
%       1, x, y, x^2, xy, y^2
%   i.e. TERMS = [0 0; 1 0; 0 1; 2 0; 1 1; 0 2].
%
%   See also adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.mewlsSolver.

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
