"""
Routines for performing the discrete hankel transform, transcribed from the MATLAB code HankelTransform
(https://documents.epfl.ch/users/l/le/leuteneg/www/MATLABToolbox/HankelTransform.html)

** NOTE: This is not finished at all at this point... **
"""
# def discrete_ht(h,R,N=None,n=None):
#    if N is None:
#       N=256
#
#    if n is None:
#       n=0
#
#    if len(n) > 1:
#       K=N
#       I=n
#    else:
#       if ~isempty(h) & isa(h,'numeric')
#          error('Need a function h(r) without kernel.');
#       end
#
#       load('dht.mat');                 % Bessel Jn rooths
#       C=c(1+n,1+N);
#       c=c(1+n,1:N);
#       r=R/C*c(:);
#       k=c(:)/R;
#       I=abs(besselj(1+n,c));
#       if n > 0
#          I(1)=1/N;                     % avoid zero - thanks to Nicolas Grisouard
#       end
#    %    I(~I)=1/N;                     % or this, but there should be only one
#       K=2*pi*R/C*I(:);
#       R=I(:)/R;
#       I=sqrt(2/C)./I;
#       I=I(:)*I.*besselj(n,c(:)/C*c);
#    end
#    if isempty(h)
#       H=h;
#    else
#       if ~isa(h,'numeric')
#          h=feval(h,r);
#       end
#       H=I*(h./R).*K;
#    end
