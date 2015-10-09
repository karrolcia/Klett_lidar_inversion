% This function is calculating extinction according to the Klett equation
% Input: 
% z  : range in km  (vector)
% P  : Lidar power (arbitary units) (vector)
% sd_P : standard deviation of P (arbitary units) (vector)
% zn1  : Lower limit of normalization range (scalar)
% zn2  : Upper limit of normalization range (zo = (zn1+zn2)/2 (scaler)
% alpha_o: Average extinction to corresponding to normalization (scaler)
% interval
% atten_beta : attenuated backscatter (Range Corrected Signal)
% 
% Output:
% 
% alpha   : Extinction coeffcient in 1/km (vector)
% d_alpha : Extimated error due to signal noise (vector)
%
% September 2013
% K. Sarna (k.sarna@tudelft.nl)

function [alphap, d_alphap] = Klett_inversion(z, Pp, sd_Pp, zn1, zn2, alpha_o)
	
%  Setting up needed arrays	  
h_res = z(2)-z(1) ;            % Range resolution in m
n = max(Pp > 0) ;              % size of i/o arrays - number of heights

% Using Range Corrected Signal RCS=P*z^2 (no ln because of the later
% change in Klett equation (without ln not necessary to calculate exp)
 
ik = find((z <= zn2) & (z >= zn1));
% ii_k = ind2sub(size(ii_k),ii_k) ;
% [i_k, j_k] = ind2sub(siz, ii_k) ;
Pp_o = mean(Pp(ik)) ;

% Find the index corresponding to zo
zo = (zn1+zn2)/2 ;
io  = find (z <= zo,1, 'last') ;
% [io, jo] = max(find((z <= zo))) ; % make sure this works well!
%Integrating RCS over the interval
int_term = zeros(1,n);
int_term(1) = 0.5 * Pp(1)*z(1)^2;

for i=2:n 
    int_term(i) = int_term(i-1) + 0.5 * (Pp(i-1) + Pp(i)) * h_res ;
end
  
int_term = (int_term-int_term(io)) ./ Pp_o ;
%int_term = trapz(z,Pp);

%Calculating extinction 

alphap = (Pp./Pp_o) ./ ((1.0/alpha_o) - 2.0 .* int_term) ;

% Calculate the Error in alpha due to the error in P

d_alphap = zeros(1,n);

 for i= 1:n    
   if (i <= io) 
    d_alphap(i) = ((alphap(i)/Pp(i)) + 2.0*alphap(i)^2 / Pp(i)*z(i)^2* h_res)^2 ...
                * sd_Pp(i)^2 ;
   elseif ((i+1) <= io)
    d_alphap(i) = d_alphap(i) + (2.0*alphap(i)^2 /Pp(i)/z(i)^2)^2 * h_res^2 ...
                * sum((z(i+1:io)).^2 * sd_Pp(i+1:io).^2) ;
%    elseif ((io+1) <= i) 
%     d_alphap(i) = d_alphap(i) + (2.0 * (alphap(i) ^2  / Pp(i) / z(i)^2))^2 ...
%                 * h_res^2 * sum((z(io+1:i).^2 * sd_Pp(io+1:i).^2)) ;
   else
    d_alphap(i) = ((alphap(i) / Pp(i)) + 2.0 *(alphap(i)^2 / Pp(i) * z(i)^2)^2) ...
                * h_res^2 * sd_Pp(i)^2 ;
   end
 end

% Adding the error due to the error in Po

sd_P_o = var(Pp(ik)) / sqrt(numel(ik)) ;  
% d_alpha = d_alpha + (alpha^2/P_o/z^2*(z[io]/alpha[io]+2*z[io]^2*res))^2*sd_P_o^2
d_alphap = d_alphap + (alphap.^2 / Pp_o / z.^2 * (z(io)/alphap(io) + 2.0*z(io)^2 * h_res)).^2 * sd_P_o^2 ;
d_alphap = sqrt(d_alphap);
