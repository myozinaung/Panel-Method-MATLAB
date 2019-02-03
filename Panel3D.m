% ************************************************************************ %
%   The 3D-Potential Flow xadded mass computation full                     %
%  	Yun Sok LEE                                                        %
% ************************************************************************ % 
clear;
load('Containership.mat');

IT = 51; nit = 50; mit = 50; nm2 = nit*mit;
nm2p1 = nm2+1;
LL = 10; MM=10; LLMM=100;

% Initialize
A    = zeros(nm2,nm2);
B    = zeros(nm2);
C    = zeros(nm2);
IPIV = zeros(nm2);

x  = zeros(IT,IT);
y  = zeros(IT,IT);
z  = zeros(IT,IT);
xp = zeros(nm2);
yp = zeros(nm2);
zp = zeros(nm2);

area = zeros(nm2);
anx  = zeros(nm2);
any  = zeros(nm2);
anz  = zeros(nm2);

u   = zeros(nm2);
v   = zeros(nm2);
w   = zeros(nm2);
pb  = zeros(nm2);
phi = zeros(nm2);
pa  = zeros(nm2);

g = zeros(nm2p1,nm2p1);

%% Ask INPUT file and prepare OUTPUT file

% 	write(*,*) 'output file name'
% 	read(*,'(a)') file2     
%       open(4,file=file2)
% %     integrated force output file name
%       open(10,file='HS3Dforce.dat',status='unknown')
% %     source distribution will be stored for field point computation at diffrent time
%       open(20,file='HS3D_source.dat',status='unknown')
%       open(4,file='field.dat')

% n =  no. of panel in longitudinal direction
% m =  no. of panel in transverse (only half-breadth)
% Hull = (x,y,z) points 

nm   = n*m;
np1  = n+1;   % node points = no. of panel + 1 (longitudinal)
mp1  = m+1;   % node points = no. of panel + 1 (transverse)
nmp1 = nm+1;
      
pi   = 3.1415927;
u0   = 1;
dudt = 1;

for i = 1:np1
    for j = 1:mp1 
        x(i,j) = Hull(((i-1)*mp1)+j,1);
        y(i,j) = Hull(((i-1)*mp1)+j,2);
        z(i,j) = Hull(((i-1)*mp1)+j,3);
    end
end
% INPUT reading finished

%% Calculation Starts
nn = 0;
for i = 1:n
    for j = 1:m
        nn = nn+1;
        
        xp(nn) = 0.25 * (x(i,j) + x(i+1,j) + x(i,j+1) + x(i+1,j+1));
        yp(nn) = 0.25 * (y(i,j) + y(i+1,j) + y(i,j+1) + y(i+1,j+1));
        zp(nn) = 0.25 * (z(i,j) + z(i+1,j) + z(i,j+1) + z(i+1,j+1));
        
        x1 = 0.5 * ((x(i+1,j) + x(i+1,j+1)) - (x(i,j) + x(i,j+1)));
        y1 = 0.5 * ((y(i+1,j) + y(i+1,j+1)) - (y(i,j) + y(i,j+1)));
        z1 = 0.5 * ((z(i+1,j) + z(i+1,j+1)) - (z(i,j) + z(i,j+1)));
        
        x2 = 0.5 * ((x(i,j) + x(i+1,j)) - (x(i,j+1) + x(i+1,j+1)));
        y2 = 0.5 * ((y(i,j) + y(i+1,j)) - (y(i,j+1) + y(i+1,j+1)));
        z2 = 0.5 * ((z(i,j) + z(i+1,j)) - (z(i,j+1) + z(i+1,j+1)));
        
        va = y1*z2 - y2*z1;
        vb = x2*z1 - x1*z2;
        vc = x1*y2 - y1*x2;
        
        area(nn) = sqrt(va^2 + vb^2 + vc^2);
        
        anx(nn) = va/area(nn);
        any(nn) = vb/area(nn);
        anz(nn) = vc/area(nn);

    end
end

for j = 1:nm
    disp([num2str(j),'th panel Boundary Condition for matrix is computed'])
    for i = 1:nm
        [uij,vij,wij] = VELO(i,j,m,x,y,z,xp,yp,zp,anx,any,anz);
        A(i,j) = uij*anx(j) + vij*any(j) + wij*anz(j);
    end
    B(j)=(-1)*anx(j);
end
 

for j = 1:nm
    for i = 1:nm
        g(j,i) = A(i,j);
    end
    g(j,nmp1) = B(j);
end
  
disp('Before matrix solver');
[g, ill] = SWEEPN(g,nm2p1,nm2p1,nm,nmp1);
disp('After matrix solver')

c =  zeros(nm,1);
HS3D_source = zeros(nm,4);
for j = 1:nm
    c(j) = g(j,nmp1);
    HS3D_source(j,1) = xp(j);
    HS3D_source(j,2) = yp(j);
    HS3D_source(j,3) = zp(j);
    HS3D_source(j,4) = c(j);
end
   
FA = 0;
FB = 0;
FT = 0;

% Header of OUTPUT file
% write(4,*) 'TITLE = "Example: Hess Smith 3D"'
% write(4,*) 'VARIABLES = "X", "Y","Z","u","v","w","ps","pu"'
% write(4,202) m,n
% 202 format(1x,'ZONE T="BIG ZONE", I=',i5,', J=',i5,', F=POINT')

OUTPUT = zeros(nm,8);
for j = 1:nm
    disp([num2str(j), ' th panel velocity and potential are computed']);
    
    u(j)   = 1;
    v(j)   = 0;
    w(j)   = 0;
    phi(j) = 0;
    
    for i = 1:nm
        [uij,vij,wij] = VELO(i,j,m,x,y,z,xp,yp,zp,anx,any,anz);
        phi0 = potent(i,j,m,x,y,z,xp,yp,zp);
        
        u(j)   = u(j)   + c(i)*uij;
        v(j)   = v(j)   + c(i)*vij;
        w(j)   = w(j)   + c(i)*wij;
        phi(j) = phi(j) + c(i)*phi0;
    end
    UVW   = sqrt(u(j)^2 + v(j)^2 + w(j)^2);
    
    pa(j) = -phi(j)*dudt; 
    pb(j) = 0.5*u0^2*(1-(UVW^2));
    
    FA = FA + (-anx(j)*pa(j)*area(j));
    FB = FB + (-anx(j)*pb(j)*area(j));
    FT = FT + ((pa(j)+pb(j))*area(j)*(-anx(j)));
    
%     write(4,1001) xp(j),yp(j),zp(j),u(j),v(j),w(j),pb(j),pa(j);
    OUTPUT(j,1) = xp(j);
    OUTPUT(j,2) = yp(j);
    OUTPUT(j,3) = zp(j);
    OUTPUT(j,4) = u(j);
    OUTPUT(j,5) = v(j);
    OUTPUT(j,6) = w(j);
    OUTPUT(j,7) = pb(j);
    OUTPUT(j,8) = pa(j);
end
 
HS3Dforce = [FA,FB,FT];
%  1001   FORMAT(8F10.5)

%     for good drawing, we try to get the original grid point data by int and ext
      
%      field point
for jjj = 1:11
    for iii = 1:21 
        yf = real(iii-1)*0.02;
        zf = real(jjj-1)*0.02;
        xf = 1.0;
        
        uf   = 1;
        vf   = 0;
        wf   = 0;
        phif = 0;
        
        for i = 1:nm
            [uij,vij,wij] = FVELO(i,iii,xf,yf,zf,m,x,y,z,anx,any,anz); % input 'iii' is added by Myo
            %   call potent(i,j,phi0,m,n);
            uf = uf + c(i)*uij;
            vf = vf + c(i)*vij;
            wf = wf + c(i)*wij;
            %	PHIf=phif+c(i)*phi0;
        end
        uvw = sqrt(uf^2 + vf^2 + wf^2);
        paf = -phif*dudt;
        pbf = 0.5*u0^2*(1-(UVW^2));
        
%         write(4,1001) xf,yf,zf,uf,vf,wf,pbf,paf
    end
end

