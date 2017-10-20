% Tahmid Mehdi
% contains the coarse-grained nonlinear system of differential equations for modelling the Master-Slave Latch
% April 2016

function dy = mslatch(t,y)
    global kA2; global kKA; global lambdap; global kB2; global kJB; global kJ2; global kK2;
    global v; global Khomodimers; global Kheterodimers; global lambda;
    syms a2 ka k2 b2 jb j2;
    AK = vpasolve([kA2(1)*(y(1,1)-2*a2-ka)^2-kA2(2)*a2-lambdap*a2==0,...
        kKA(1)*(y(4,1)-2*k2-ka)*(y(1,1)-2*a2-ka)-kKA(2)*ka-lambdap*ka==0,...
        kK2(1)*(y(4,1)-2*k2-ka)^2-kK2(2)*k2-lambdap*k2==0],...
        [a2,ka,k2], [0 Inf; 0 Inf; 0 Inf]);
    A2=double(max(AK.a2)); KA=double(max(AK.ka)); K2=double(max(AK.k2));
    BJ = vpasolve([kB2(1)*(y(2,1)-2*b2-jb)^2-kB2(2)*b2-lambdap*b2==0,...
        kJB(1)*(y(3,1)-2*j2-jb)*(y(2,1)-2*b2-jb)-kJB(2)*jb-lambdap*jb==0,...
        kJ2(1)*(y(3,1)-2*j2-jb)^2-kJ2(2)*j2-lambdap*j2==0],...
        [b2,jb,j2], [0 Inf; 0 Inf; 0 Inf]);
    B2=double(max(BJ.b2)); JB=double(max(BJ.jb)); J2=double(max(BJ.j2));
    OA2 = (1+A2/Khomodimers(1))^(-2);
    OB2 = (1+B2/Khomodimers(2))^(-2);
    OJ2 = (1+J2/Khomodimers(3))^(-2);
    OK2 = (1+K2/Khomodimers(4))^(-2);
    %OKA = (1+KA/(Kheterodimers(1)*(1+JB/Kheterodimers(2))))^(-1);
    %OJB = (1+JB/(Kheterodimers(2)*(1+KA/Kheterodimers(1))))^(-1);
    OKA = (1+KA/Kheterodimers(1))^(-1);
    OJB = (1+JB/Kheterodimers(2))^(-1);
    OA2act = (A2/Khomodimers(1)/(1+A2/Khomodimers(1)))^2;
    dy(1,1) = v(1)*OB2*OKA-lambda*y(1,1);
    dy(2,1) = v(2)*OA2*OJB-lambda*y(2,1);
    dy(3,1) = v(3)*OK2*OKA*OA2-lambda*y(3,1);
    dy(4,1) = v(4)*OJ2*OJB-lambda*y(4,1);
    % flip the switch
    if t > 60
        v(3) = 11.6;
        v(4) = 0;
    end
end