function [A, ILL] = SWEEPN(A,ND,MD,N,M)
%  SUBROUTINE OF SWEEP-OUT METHOD    ( BY GAUSSIAN ELIMINATION METHOD )
% real*8 A(ND,MD)

if(N > ND || M > MD || N >= M || N <= 1) % goto 200
    ILL = - 90000;
    disp(['Input Data Error: ',ND,' ',MD,' ',N,' ',M]);
    return;
end

ILL  = 0;
EPS  = 1.0;
PMAX = 0;
%% ++++++  NORMARIZATION  ++++++
for I = 1:N
    Z = 0.0;
    for J = 1:N
        W = abs(A(I,J));
        if (W > Z) 
        Z = W;
        end
    end
    if(Z < 1.E-40) % goto 100
        ILL = I;
        disp(['Matrix is singular: ',ILL]);
        return;
    end
    W = 1.0/Z;
    for J = 1:M
        A(I,J) = A(I,J)*W;
    end
end
    
%%      ++++++  PARTIAL  PIVOTING   ( ROW  EXCHANGE )  ++++++
N1 = N-1;

for I=1:N1
    PIV = 0.0;

    for K = I:N
        W = abs(A(K,I));
        if (W < PIV)
            continue;
        end
        PIV = W;
        L = K;
    end

    if (PIV < 1.E-8*PMAX) % goto 110
        ILL = I;
        disp(['Matrix is nearly singular: ',ILL]);
        return;
    end

    if (PIV < 1.E-5*PMAX) 
        ILL = ILL+1;
    end
    if (PIV < EPS) 
        EPS = PIV;
    end
    if (PIV > PMAX) 
        PMAX = PIV;
    end

    if (L == I) % go to 13
        % ++++++  FORWARD  REDUCTION  ++++++
        L = I+1; % Start of Block13
        Z = 1.0/A(I,I);
        for J = L:M
            A(I,J) = A(I,J)*Z;
        end

        for K = L:N % loop 15
            W = -A(K,I);
            if( abs(W) < 1.E-8*PMAX) 
                continue; % loop 15
            end
            for J = L:M
                A(K,J) = A(K,J) + A(I,J)*W;
            end
        end % end of Block13
        
    else
        for J = I:M % loop to skip if L = I
            W      = A(L,J);
            A(L,J) = A(I,J);
            A(I,J) = W;
        end

        % ++++++  FORWARD  REDUCTION  ++++++
        L = I+1; %  Start of Block13
        Z = 1.0/A(I,I);
        for J = L:M
            A(I,J) = A(I,J)*Z;
        end

        for K = L:N % loop 15
            W = -A(K,I);
            if( abs(W) < 1.E-8*PMAX) 
                continue; % loop 15
            end
            for J = L:M
                A(K,J) = A(K,J) + A(I,J)*W;
            end
        end  % end of Block13
    end
    
end
   %%%%%%%%%%%
   
I = N;
W = abs(A(N,N));
if (W < 1.E-8*PMAX) % goto 110
    ILL = I;
    disp(['Matrix is nearly singular: ',ILL]);
    return;
end

if (W < 1.E-5*PMAX) 
    ILL = ILL+1;
end
if (W < EPS) 
    EPS = W;
end
if (W > PMAX) 
    PMAX = W;
end

LB = N+1;
for J = LB:M
    A(N,J) = A(N,J)/A(N,N);
end

%% ++++++  BACKWARD  SUBSTITUTION  ++++++
for K = 2:N
    I  = N - K+1;
    JB = I+1;
    for J = JB:N
        W = -A(I,J);
        for L = LB:M
            A(I,L) = A(I,L) + A(J,L)*W;
        end
    end
end

if (ILL == 0) 
    return;
end
          
%% ++++++  ERROR  MESAGE  ++++++

EPS = EPS/PMAX;
disp(['The solution is inaccurate: ',ILL,' ',EPS]);
return;

% WRITE(6,99) ILL,EPS
% 99     FORMAT(/////,10X,36H***  THE SOLUTION IS INACCURATE  ***, 10X,
% &      5HILL = I4,10X,5HEPS = E10.3)
% RETURN

% % *
% 100     ILL=I
% WRITE(6,101) ILL
% 101     FORMAT(/////,10X,28H***  MATRIX IS SINGULAR  ***, 10X,
% &      6HAT THE,I4, 7H-TH ROW/////)
% RETURN

% % *
% 110   ILL=I
% WRITE(6,111) ILL
% 111   FORMAT(/////,10X,35H***  MATRIX IS NEARLY SINGULAR  ***, 10X,
% &    6HAT THE,I4,7H-TH ROW/////)
% RETURN

% % *
% 200 ILL=-90000
% WRITE(6,201) ND,MD,N,M
% 201 FORMAT(///,10X,28H***  INPUT  DATA  ERROR  ***, 20X,11HARGUMENTS = 4(I5,1H,)/////)
% RETURN

end % end of function
