function  [CH, JCH] = gsp_jackson_cheby_coeff(a, b, lambda_range, m)
% Compute the coefficients of the polynomial approximation of an ideal
% bandpass filter.
%
%   [CH, JCH] = jackson_cheby_poly_coefficients(a, b, lambda_range, m)
%
% Ouputs:
%   - CH contains the list of chebyshev coefficients
%   - JCH contains the list of Jackson chebyshev coefficients
%
% Inputs:
%   - a and b indicates the start and end of the bandpass interval.
%   - lambda_range is typically [0 G.lmax] in graph signal processing and
%   indicates the full interval on which the filter is defined.
%   - m is the order of the polynomial approximation.
%
%
% This code is implemented using the methods described in:
% [1] E. D. Napoli, E. Polizzo, and Y. Saad, "Efficient estimation of 
% eigenvalue counts in an interval," arXiv:1308.4275, 2013.
% [2] D. K. Hammond, P. Vandergheynst, and R. Gribonval, "Wavelets on 
% graphs via spectral graph theory,” Appl. Comput. Harmon. Anal., vol. 30, 
% no. 2, pp. 129–150, 2011
%
% The chebyshev coefficients can also be computed using the function 
% gsp_cheby_coeff in the GSP toolbox. Details about the GSP toolbox can be
% found in:
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
%   G. Puy, N. Tremblay, R. Gribonval and P. Vandergheynst. "Random 
%   sampling of bandlimited signals on graphs", ACHA, 2016.


% scaling and translation coefficients compared to the classical interval
% of Chebychev polynomials [-1,1] :
a1 = (lambda_range(2)-lambda_range(1))/2;
a2 = (lambda_range(1)+lambda_range(2))/2;

% scale the boundaries of the band pass according to lrange:
a=(a-a2)/a1;
b=(b-a2)/a1;

% compute Cheby coef:
CH(1)=(1/pi)*(acos(a)-acos(b));
for j=2:m+1
    CH(j)=(2/(pi*(j-1)))*(sin((j-1)*acos(a))-sin((j-1)*acos(b)));
end

% compute Jackson coef: specific to ideal bandpass functions
alpha=pi/(m+2);
for j=1:m+1
    gamma_JACK(j)=(1/sin(alpha))*((1-(j-1)/(m+2))*sin(alpha)*cos((j-1)*alpha)+(1/(m+2))*cos(alpha)*sin((j-1)*alpha));
end

% compute Jackson-Cheby coef:
JCH=CH.*gamma_JACK;

% to be in adequation with gsp_cheby_op.m :
JCH(1)=JCH(1)*2;
CH(1)=CH(1)*2;

JCH=JCH';
CH=CH';