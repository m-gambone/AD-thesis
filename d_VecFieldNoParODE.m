% Generated by ADiMat 0.6.2-5299 (cf36599d)
% © 2001-2008 Andre Vehreschild <vehreschild@sc.rwth-aachen.de>
% © 2009-2015 Johannes Willkomm <johannes@johannes-willkomm.de>
% TU Darmstadt, 64289 Darmstadt, Germany
% Visit us on the web at http://www.adimat.de/
% Report bugs to johannes@johannes-willkomm.de
%
%                             DISCLAIMER
% 
% ADiMat was prepared as part of an employment at the Institute for Scientific Computing,
% RWTH Aachen University, Germany and at the Institute for Scientific Computing,
% TU Darmstadt, Germany and is provided AS IS. 
% NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL REPUBLIC OF GERMANY
% NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY, NOT THE TU DARMSTADT,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
% OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
% OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE
% WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
%
% Parameters:
%  - dependents=der_S, der_I, der_E
%  - independents=S, I, E
%  - inputEncoding=ISO-8859-1
%
% Functions in this file: d_VecFieldNoParODE
%

function [d_der_S der_S d_der_I der_I d_der_E der_E] = d_VecFieldNoParODE(d_S, S, d_I, I, d_E, E, tau)
%Funzione per il calcolo del campo vettoriale associato al modello SIR con vaccinazione
   n = length(E);
   mu = 1 / (75 * 365);
   p0 = 0.75;
   nu = 1 / 7;
   beta = 20 * (mu + nu);
   a = n / tau;
   c = 300;
   tmpda2 = 1 - p0;
   d_tmpca1 = adimat_opdiff_mult_left(c, adimat_opdiff_subsref(d_E, struct('type', '()', 'subs', {{n}})), E(n));
   tmpca1 = c * E(n);
   [d_p1 p1] = adimat_diff_min2(d_tmpca1, tmpca1, d_zeros(tmpda2), tmpda2);
   k = 1;
   d_g = adimat_opdiff_mult_left(k, d_I, I);
   g = k * I;
   d_tmpca5 = adimat_opdiff_mult_left(beta, d_S, S);
   tmpca5 = beta * S;
   d_tmpca4 = adimat_opdiff_mult(d_tmpca5, tmpca5, d_I, I);
   tmpca4 = tmpca5 * I;
   d_tmpca3 = adimat_opdiff_mult_left(mu, d_S, S);
   tmpca3 = mu * S;
   d_tmpca2 = adimat_opdiff_sum(-d_p1, d_zeros(1 - p0));
   tmpca2 = 1 - p0 - p1;
   d_tmpca1 = adimat_opdiff_mult_left(mu, d_tmpca2, tmpca2);
   tmpca1 = mu * tmpca2;
   d_der_S = adimat_opdiff_sum(d_tmpca1, -d_tmpca3, -d_tmpca4);
   der_S = tmpca1 - tmpca3 - tmpca4;
   d_tmpca2 = adimat_opdiff_mult_left(beta, d_S, S);
   tmpca2 = beta * S;
   d_tmpca1 = adimat_opdiff_sum(d_tmpca2, d_zeros(-mu - nu));
   tmpca1 = tmpca2 - mu - nu;
   d_der_I = adimat_opdiff_mult(d_I, I, d_tmpca1, tmpca1);
   der_I = I * tmpca1;
   d_tmpca2 = adimat_opdiff_mult_left(a, adimat_opdiff_subsref(d_E, struct('type', '()', 'subs', {{1}})), E(1));
   tmpca2 = a * E(1);
   d_tmpca1 = adimat_opdiff_mult_left(a, d_g, g);
   tmpca1 = a * g;
   d_der_E = adimat_opdiff_subsasgn([], struct('type', {'()'}, 'subs', {{1}}), adimat_opdiff_sum(d_tmpca1, -d_tmpca2));
   der_E(1) = tmpca1 - tmpca2;
   for i=2 : n
      d_tmpca3 = adimat_opdiff_mult_left(a, adimat_opdiff_subsref(d_E, struct('type', '()', 'subs', {{i}})), E(i));
      tmpca3 = a * E(i);
      tmpda2 = i - 1;
      d_tmpca1 = adimat_opdiff_mult_left(a, adimat_opdiff_subsref(d_E, struct('type', '()', 'subs', {{tmpda2}})), E(tmpda2));
      tmpca1 = a * E(tmpda2);
      d_der_E = adimat_opdiff_subsasgn(d_der_E, struct('type', {'()'}, 'subs', {{i}}), adimat_opdiff_sum(d_tmpca1, -d_tmpca3));
      der_E(i) = tmpca1 - tmpca3;
   end
end
