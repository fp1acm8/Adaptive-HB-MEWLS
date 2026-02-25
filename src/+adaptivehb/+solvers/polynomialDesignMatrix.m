function [A, terms] = polynomialDesignMatrix(data, degree)
%POLYNOMIALDESIGNMATRIX Build a 2D polynomial design matrix.
%   [A, TERMS] = POLYNOMIALDESIGNMATRIX(DATA, DEGREE) returns the design
%   matrix A for all bivariate monomials x^px * y^py whose total degree
%   px + py does not exceed DEGREE.
%
%   In the notation of Brugnano et al. (2024), this constructs the matrix
%   Phi for the bivariate polynomial space Pi_d^2.  Each column of Phi
%   corresponds to one monomial basis function phi_j(x,y) = x^p * y^q.
%
%   Inputs:
%       DATA   - N-by-2 matrix of coordinates (typically in [0,1]^2).
%       DEGREE - non-negative integer setting the maximum total degree.
%
%   Outputs:
%       A     - N-by-M design matrix where M = (degree+1)*(degree+2)/2.
%               Each column corresponds to one monomial term.
%       TERMS - M-by-2 matrix of exponent pairs [px, py].
%
%   The terms are enumerated in graded lexicographic order.  For example,
%   with degree = 2 the columns of A correspond to:
%       1, x, y, x^2, xy, y^2
%   i.e. TERMS = [0 0; 1 0; 0 1; 2 0; 1 1; 0 2].
%
%   See also adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.mewlsSolver.

% --- Input validation (MATLAB R2020b+ arguments block) ---------------
% data must have exactly 2 columns; degree must be a non-negative integer.
arguments
    data (:, 2) double
    degree (1, 1) double {mustBeInteger, mustBeNonnegative}
end

% Number of data points (rows of the design matrix).
nPts = size(data, 1);

% Generate all exponent pairs (p, q) with p+q <= degree.
% The total number of terms is M = (d+1)(d+2)/2 for degree d.
terms = generate_terms(degree);

% Pre-allocate the design matrix Phi with ones. The constant term (p=q=0)
% column is already correct; the other columns will be overwritten below.
A = ones(nPts, size(terms, 1));

% Extract the x and y coordinate vectors for readability.
x = data(:, 1);  % N-by-1, first spatial coordinate
y = data(:, 2);  % N-by-1, second spatial coordinate

% Fill each column of Phi with the monomial x^p * y^q.
for i = 1:size(terms, 1)
    px = terms(i, 1);  % exponent of x for this term
    py = terms(i, 2);  % exponent of y for this term

    if px == 0 && py == 0
        % Constant term phi = 1 (already set during pre-allocation).
        A(:, i) = 1;
    else
        % General monomial: phi_i(x,y) = x^px * y^py.
        % Element-wise power and multiplication across all N points.
        A(:, i) = (x .^ px) .* (y .^ py);
    end
end

end


% =====================================================================
%  Local function: generate_terms
%  Enumerates all bivariate exponent pairs (px, py) with px+py <= d.
% =====================================================================
function terms = generate_terms(degree)
%GENERATE_TERMS List all monomial exponents in graded lex order.
%   Returns an M-by-2 matrix where M = (d+1)(d+2)/2.
%   Ordering: for each total degree k = 0,1,...,d, we enumerate px from
%   k down to 0 (equivalently px from 0..d, py from 0..d-px).

terms = [];  % will grow to M-by-2

% Outer loop over x-exponent, inner over y-exponent.
% For each px, py ranges from 0 to (degree - px) so that px+py <= d.
for px = 0:degree
    for py = 0:(degree - px)
        terms(end + 1, :) = [px, py]; %#ok<AGROW>
    end
end
% Result: terms is M-by-2 with M = (degree+1)*(degree+2)/2.
end
