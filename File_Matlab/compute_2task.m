clear

syms teta d a alfa real;
syms teta1 teta2 teta4 teta3 teta5 teta6 teta7 real;
a1=0; a2=0;
a3=0; a4=0; a5=0; a6=0; a7=0;
alfa1= -pi/2; 
alfa2= pi/2; 
alfa3= pi/2; 
alfa4= -pi/2;
alfa5= -pi/2;
alfa6= pi/2;
alfa7= 0;
d1=0.34; 
d2=0; 
d3=0.4;
d4=0;
d5=0.4;
d6=0;
d7=0.86;
%d7=0.126;
%k0=10;



A=[cos(teta) -sin(teta)*cos(alfa) sin(teta)*sin(alfa) a*cos(teta); 
    sin(teta) cos(teta)*cos(alfa) -cos(teta)*sin(alfa) a*sin(teta);
    0 sin(alfa) cos(alfa) d;
    0 0 0 1];

A01=subs(A, [alfa teta a d], [alfa1 teta1 a1 d1] );
A01=simplify(A01);
A12=subs(A, [alfa teta a d], [alfa2 teta2 a2 d2] );
A12=simplify(A12);
A23=subs(A, [alfa teta a d], [alfa3 teta3 a3 d3] );
A23=simplify(A23);
A34=subs(A, [alfa teta a d], [alfa4 teta4 a4 d4] ); 
A34=simplify(A34);
A45=subs(A, [alfa teta a d], [alfa5 teta5 a5 d5] );
A45=simplify(A45);
A56=subs(A, [alfa teta a d], [alfa6 teta6 a6 d6] );
A56=simplify(A56);
A67=subs(A, [alfa teta a d], [alfa7 teta7 a7 d7] );
A67=simplify(A67);

Tbe=A01*A12*A23*A34*A45*A56*A67;
Tbe=simplify(Tbe);

p_q = Tbe(1:3, 4);

obj=[0.75 0.75 0.2]';
p = p_q - obj;
p=simplify(p);
dist = norm(p);
dist=simplify(dist);
wdq=jacobian(dist, [teta1, teta2, teta3, teta4,teta5, teta6, teta7]);
wdq=simplify(wdq);
wdq = wdq';
