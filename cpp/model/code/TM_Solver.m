function [r J] = TM_Solver()

    nbpop = 2 ;
    Iext = [ 1 .5 ] ;
    J = [ .44 -0.28 ; 0.67 -0.34 ] ;

% Iext = [0.1 0.1];
    % J = [1.92 -0.18; 0.25 -0.20] ;

    U = [ .3 1 ; 1 .3 ] ;
    Tfac = [ 1 0 ; 0 .03 ] ;
    Trec = [ .03 0 ; 0 1 ] ;

    r0 = rand(1,nbpop) ;
    r = -r0 ;
    options = optimset('Display','off') ;
    while(any(r<0))
        r0 = rand(1,nbpop) ;
        J = rand(nbpop,nbpop) ;
        J(:,2) = -J(:,2) ;
        r = fsolve(@SelfCstEq,r0,options) ;
    end

    disp(sprintf('rates'))
      disp(sprintf('%.3f ', r))
      disp(sprintf('error'))
      disp(sprintf( '%.3f ', SelfCstEq(r)) )
     
    function Eq = SelfCstEq(r)
        Eq = [] ;
        for i=1:nbpop
            Eq_i = Iext(i) ;
            for j=1:nbpop
                
                u_ij = U(i,j) * ( 1 + Tfac(i,j) * r(i) ) / ( 1 + U(i,j) * Tfac(i,j) * r(i) ) ; 
                x_ij = 1 / ( 1 + Trec(i,j) * u_ij * r(i) ) ; 
                
                Eq_i = Eq_i + J(i,j) * u_ij * x_ij * r(j) ;
            end
            Eq = [Eq ; Eq_i] ;
        end
    end
    
end
