% Wigner3j.m by David Terr, Raytheon, 6-17-04
% Compute the Wigner 3j symbol using the Racah formula [1]. 
% updated by Adi Rangan 20181225; changing factorial to (safer) dfactorial and lfactorial. ;
% test with: ;
%{
  n_iteration = 1024; e_sum = 0;
  for niteration = 1:n_iteration;
  j1 = round(15*rand());
  j2 = round(15*rand());
  j3 = abs(j1-j2) + floor(((j1+j2)-abs(j1-j2))*rand());
  m1 = j1*round(2*rand()-1);
  m2 = j2*round(2*rand()-1);
  m3 = j3*round(2*rand()-1);
  e_sum = e_sum + abs(wigner3j_0(j1,j2,j3,m1,m2,m3) - wigner3j_1(j1,j2,j3,m1,m2,m3));
  end;%for niteration = 1:n_iteration;
  disp(sprintf(' %% n_iteration %d e_sum: %0.16f',n_iteration,e_sum));
  %}
function wigner = Wigner3j_1(j1,j2,j3,m1,m2,m3)
% error checking
if ( 2*j1 ~= floor(2*j1) || 2*j2 ~= floor(2*j2) || 2*j3 ~= floor(2*j3) ...
        || 2*m1 ~= floor(2*m1) || 2*m2 ~= floor(2*m2) || 2*m3 ~= floor(2*m3) )
    %error('All arguments must be integers or half-integers.');
    wigner = 0.0d0;
    return;
end
if ( j1 - m1 ~= floor ( j1 - m1 ) )
    %error('2*j1 and 2*m1 must have the same parity');
    wigner = 0.0d0;
    return;
end
if ( j2 - m2 ~= floor ( j2 - m2 ) )
    %error('2*j2 and 2*m2 must have the same parity');
    wigner = 0.0d0;
    return;
end
if ( j3 - m3 ~= floor ( j3 - m3 ) )
    %error('2*j3 and 2*m3 must have the same parity');
    wigner = 0.0d0;
    return;
end
if j3 > j1 + j2 || j3 < abs(j1 - j2)
    %error('j3 is out of bounds.');
    wigner = 0.0d0;
    return;
end
if abs(m1) > j1
    %error('m1 is out of bounds.');
    wigner = 0.0d0;
    return;
end
if abs(m2) > j2
    %error('m2 is out of bounds.');
    wigner = 0.0d0;
    return;
end
if abs(m3) > j3
    %error('m3 is out of bounds.');
    wigner = 0.0d0;
    return;
end
if (m1+m2+m3~=0)
    wigner = 0.0d0;
    return;
end
t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;
tmin = max( 0, max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );
wigner = 0;
for t = tmin:tmax
    wigner = wigner + (-1)^t / exp( lfactorial(t) + lfactorial(t-t1) + lfactorial(t-t2) ...
        + lfactorial(t3-t) + lfactorial(t4-t) + lfactorial(t5-t) );
end
wigner = wigner * (-1)^(j1-j2-m3) ...
  * sqrt( exp( lfactorial(j1+j2-j3) + lfactorial(j1-j2+j3) + lfactorial(-j1+j2+j3) - lfactorial(j1+j2+j3+1)...
	       + lfactorial(j1+m1) + lfactorial(j1-m1) + lfactorial(j2+m2) + lfactorial(j2-m2) + lfactorial(j3+m3) + lfactorial(j3-m3) ) );

