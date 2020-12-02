function redsys = reduceSys(sys,V,W)
    J = W'*sys.J*W;
    R = W'*sys.R*W;
    Q = V'*sys.Q*V;
    [J_out,R_out] = correctSymmetry(J,R);
    B = W'*sys.B;
    
    redsys = phs(J_out,R_out,Q,B);