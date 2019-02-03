function phi0 = potent(iii,j1,m,x,y,z,xp,yp,zp)

LL = 10; MM=10; LLMM=100;

xx = zeros(LL+1,MM+1);
yy = zeros(LL+1,MM+1);
zz = zeros(LL+1,MM+1);

xxp = zeros(LLMM);
yyp = zeros(LLMM);
zzp = zeros(LLMM);
ALM = zeros(LLMM);
  
% k = INT(iii/m) block
k = iii/m;
if isreal(k)
    if abs(k) < 1
        k = 0;
    else
        k = sign(k)*floor(abs(k));
    end
else
    k  = iii/m;
end % end of block

k2 = iii - m*k;

if (k2 == 0)
    i = k;
    j = m;
else
    i = k+1;
    j = k2;
end

for ii = 1:LL+1
    xx(ii,1) = x(i,j) + real(ii-1)*(x(i+1,j)-x(i,j))/real(LL);
    yy(ii,1) = y(i,j) + real(ii-1)*(y(i+1,j)-y(i,j))/real(LL);
    zz(ii,1) = z(i,j) + real(ii-1)*(z(i+1,j)-z(i,j))/real(LL);
    
    xx(ii,MM+1) = x(i,j+1) + real(ii-1)*(x(i+1,j+1)-x(i,j+1))/real(LL);
    yy(ii,MM+1) = y(i,j+1) + real(ii-1)*(y(i+1,j+1)-y(i,j+1))/real(LL);
    zz(ii,MM+1) = z(i,j+1) + real(ii-1)*(z(i+1,j+1)-z(i,j+1))/real(LL);
end

for ii = 1:LL+1
    for jj = 1:MM+1
        xx(ii,jj) = xx(ii,1) + real(jj-1)*(xx(ii,MM+1)-xx(ii,1))/real(MM);
        yy(ii,jj) = yy(ii,1) + real(jj-1)*(yy(ii,MM+1)-yy(ii,1))/real(MM);
        zz(ii,jj) = zz(ii,1) + real(jj-1)*(zz(ii,MM+1)-zz(ii,1))/real(MM);
    end
end

nnn = 0;
for ii = 1:LL
    for jj = 1:MM  	
        nnn = nnn+1;
        
        xxp(nnn) = 0.25*(xx(ii,jj) + xx(ii+1,jj) + xx(ii,jj+1) + xx(ii+1,jj+1));
        yyp(nnn) = 0.25*(yy(ii,jj) + yy(ii+1,jj) + yy(ii,jj+1) + yy(ii+1,jj+1));
        zzp(nnn) = 0.25*(zz(ii,jj) + zz(ii+1,jj) + zz(ii,jj+1) + zz(ii+1,jj+1));
        
        xx1 = 0.5*((xx(ii+1,jj)+ xx(ii+1,jj+1)) - (xx(ii,jj)+ xx(ii,jj+1)));
        yy1 = 0.5*((yy(ii+1,jj)+ yy(ii+1,jj+1)) - (yy(ii,jj)+ yy(ii,jj+1)));
        zz1 = 0.5*((zz(ii+1,jj)+ zz(ii+1,jj+1)) - (zz(ii,jj)+ zz(ii,jj+1)));
        
        xx2 = 0.5*((xx(ii,jj)+ xx(ii+1,jj)) - (xx(ii,jj+1)+ xx(ii+1,jj+1)));
        yy2 = 0.5*((yy(ii,jj)+ yy(ii+1,jj)) - (yy(ii,jj+1)+ yy(ii+1,jj+1)));
        zz2 = 0.5*((zz(ii,jj)+ zz(ii+1,jj)) - (zz(ii,jj+1)+ zz(ii+1,jj+1)));
        
        vva = yy1*zz2 - yy2*zz1;
        vvb = xx2*zz1 - xx1*zz2;
        vvc = xx1*yy2 - yy1*xx2;
        
        ALM(nnn) = sqrt((vva)^2 + (vvb)^2 + (vvc)^2);
    end
end

phi0 = 0;

for kk = 1:1
    phi00 = 0;
    for k = 1:LLMM
        if (kk == 1)   
            xsa = xp(j1) - xxp(k);
            ysa = yp(j1) - yyp(k);
            zsa = zp(j1) - zzp(k);
            
            RRR = sqrt(xsa^2 + ysa^2 + zsa^2);
            phi00 = phi00 + (-1)/(4*pi*RRR)*ALM(k);
            
        elseif (kk == 2)
            xsa = xp(j1) - xxp(k);
            ysa = yp(j1) + yyp(k);
            zsa = zp(j1) - zzp(k);
            
            RRR = sqrt(xsa^2 + ysa^2 + zsa^2);
            phi00 = phi00 + (-1)/(4*pi*RRR)*ALM(k);
            
        elseif(kk == 3)  
            xsa = xp(j1)-xxp(k);
            ysa = yp(j1)-yyp(k);
            zsa = zp(j1)+zzp(k);
            
            RRR = sqrt(xsa^2 + ysa^2 + zsa^2);
            phi00 = phi00 + (-1)/(4*pi*RRR)*ALM(k);
            
        else
            xsa = xp(j1) - xxp(k);
            ysa = yp(j1) + yyp(k);
            zsa = zp(j1) + zzp(k);
            
            RRR = sqrt(xsa^2 + ysa^2 + zsa^2);
            phi00 = phi00 + (-1)/(4*pi*RRR)*ALM(k);
        end
    end
    phi0 = phi00 + phi0;
end

return;
end
