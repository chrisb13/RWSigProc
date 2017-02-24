function  z=sombrero(v,h,L,T);

% function  z=sombrero(v,h,L,T);
% produces a sombrero-like shaped matrix z, used for FIR filtering
% z is v rows by h columns
% L horizontal wavelength in matrix units, may be negative, default v
% T vertical wavelength in matrix units, positive, default h

if nargin==0,
 h=5; h=v;L=v;T=v;
elseif nargin==1,
 h=v;L=v;T=v;
elseif nargin==2,
 L=h;T=v;
elseif nargin==3,
 T=v;
end

if(mod(v,2)==0),disp('Warning: v is even.');end
if(mod(h,2)==0),disp('Warning: h is even.');end

% get the gaussian envelope

s=2; % s controls the value at the borders
     % s is the sigma of the normal distribution

[x,y]=meshgrid(0.5:h-.5,0.5:v-.5);
x=(pi*(x-h/2)/(h*s)).^2; y=(pi*(y-v/2)/(v*s)).^2;
g=exp(-0.5*(x+y)) ./ (sqrt(2*pi) .* s);
g=g-min(g(:));
g=g/max(g(:));

% get the cosine surface and multiply by the gaussian

[x,y]=meshgrid( ((1:h)-(h+1)/2) , ((1:v)-(v+1)/2) );
z=g.*cos((2*pi/L)*x - (2*pi/T)*y);

% accept negative values, normalize for zero sum
sz=sum(z(:));s1=sum(ones(size(z(:))));
z=z-sz/s1;
sa=sum(abs(z(:)));
z=z/sa;