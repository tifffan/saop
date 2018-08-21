function [lambda_k, cum_coh] = gsp_estimate_lk(G, k, param)
% Function that estimate the kth smallest eigenvalue of a graph Laplacian
%
%   [lambda_k, cum_coh] = gsp_estimate_lk(G, k, param)
%
% Ouputs:
%   - lambda_k: the estimated kth smallest eigenvalue of the graph
%   Laplacian.
%   - cum_coh: the estimated *squared* graph local cumulative coherences.
%
% Inputs:
%   - G is a matlab graph structure (see GSP toolbox) representing a graph
%   - k is the index of the eigenvalue to estimate. The eigenvalues are
%   ordered from the smallest one to the largest one.
%   - param is a matlab sructure containing the options to estimate the kth
%   smallest eigenvalue. The fields are:
%       - param.nb_estimation: number of times that the complete algorithm is
%       run. The final result consists of the average of all results, to
%       improve the estimation. (Default: 1)
%       - param.nb_features: number of random vectors used to estimate 
%       lambda_k. (Default: 2*round(log(G.N))
%       - param.verbose: display information or not about the evolution of
%       the algorithm. (Default: 0 - Use 1 to display informations)
%       -  param.epsilon: tolerance to stop the algorithm. (Default: 1e-2)
%       -  param.jackson: 0 to use chebyshev polynomials to approximate the 
%       ideal filter; 1 to use Jackson chebyshev polynomials. (Default: 0)
%       - param.interval: initial interval to start the dichotomy to find 
%       lambda_k. (Default: [0, G.lmax])
%
%
% The algorithms is detailled in:
% G. Puy, N. Tremblay, R. Gribonval and P. Vandergheynst. "Random sampling
% of bandlimited signals on graphs", arxiv:1511.05118, 2016.
%
% Details about the GSP toolbox can be found in:
% N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst, 
% and D. K. Hammond, "Gspbox: A toolbox for signal processing on graphs," 
% arXiv:1408.5781, 2014
% 
% Copyright (c) 2016 G. Puy, N. Tremblay
%
% This file is part of the GraphSamplingBox
%
% The GraphSamplingBox is free software: you can redistribute it and/or 
% modify it under the terms of the GNU Affero General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%
% The GraphSamplingBox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this toolbox please kindly cite
%  G. Puy, N. Tremblay, R. Gribonval and P. Vandergheynst. "Random 
%   sampling of bandlimited signals on graphs", ACHA, 2016.

%
if ~isfield(G, 'lmax')
    G = gsp_estimate_lmax(G);
end

% Parameters
if nargin < 3
    param = struct;
end
if ~isfield(param, 'nb_estimation'),
    param.nb_estimation = 1;
end
if ~isfield(param, 'nb_features'),
    param.nb_features = 2*round(log(G.N));
end
if ~isfield(param, 'epsilon'),
    param.epsilon = 1e-2;
end
if ~isfield(param, 'verbose'),
    param.verbose = 0;
end
if ~isfield(param, 'jackson'),
    param.jackson = 0;
end
if ~isfield(param, 'interval'),
    param.interval = [0, G.lmax];
end

% Initialisation
if nargout>=2
    norm_Uk = zeros(G.N, param.nb_estimation);
end
lambda_k_est = zeros(param.nb_estimation, 1);

% Perform nb_estimation on different of set feature vectors
for ind_est = 1:param.nb_estimation
    
    % Random signals (fixed for one estimation)
    Sig = randn(G.N, param.nb_features)*1/sqrt(param.nb_features);
    
    % Search by dichotomy
    counts = 0;
    lambda_min = param.interval(1);
    lambda_max = param.interval(2);
    while counts~=k || (lambda_max - lambda_min)/lambda_max > param.epsilon
        % Middle of the interval
        lambda_mid = (lambda_min+lambda_max)/2;
        % Filter
        order = 30; % Order of the polynomial
        [ch, jch] = gsp_jackson_cheby_coeff(0, lambda_mid, ...
            [0 G.lmax], order);
        % Filtering
        if param.jackson
            X = gsp_cheby_op(G, jch(:), Sig);
        else
            X = gsp_cheby_op(G, ch(:), Sig);
        end
        % Check results
        counts = round(sum(X(:).^2));
        if counts>=k, lambda_max = lambda_mid;
        else lambda_min = lambda_mid; end
        if param.verbose,
            fprintf(' ... estimating lambda_k: %i, %e, %e\n', ...
                counts, lambda_min, lambda_max);
        end
    end
    
    % Store result
    lambda_k_est(ind_est) = (lambda_min+lambda_max)/2;
    if param.verbose,
        fprintf('%i estimation of lambda_k: %e\n', ind_est, ...
            lambda_k_est(ind_est));
    end
    if nargout>=2
        norm_Uk(:, ind_est) = sum(X.^2, 2);
    end
    
end

% Final estimation
G.lk = mean(lambda_k_est);
lambda_k = G.lk;
if param.verbose,
        fprintf('Final estimation of lambda_k: %e\n', ...
            lambda_k);
end
if nargout>=2
    cum_coh = mean(norm_Uk, 2);
end

end
