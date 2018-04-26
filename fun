clear;
pl=0;

function [x,y]=ras(z,K,n)

    disp(1);
    nv=0.5;
   //disp(nv);
    disp(2);
    fun=0;
    for i=1:n
        fun=fun+z(i)*(K(i)-1)/(nv*(K(i)-1)+1);
    end
    //disp(fun);
   
    j=0;
    //fd = mopen('C:\Users\Александр\Documents\учеба\диплом 2\text.txt','wt');
    while (abs(fun)>1e-6)// & (nv>0) & (nv<1)
        //   mfprintf(fd, '%d %d %d\n', j, fun, nv);
       // disp(j)
        j=j+1;
        //    disp(f)
        //    disp(nv)
        fun1=0;
        for i=1:n
            fun1=fun1-z(i)*(K(i)-1)^2/(nv*(K(i)-1)+1)^2
        end    
       // disp(nv);
        nv=nv-fun/fun1;
        fun=0;
         disp(7);
        for i=1:n
            fun=fun+z(i)*(K(i)-1)/(nv*(K(i)-1)+1)
        end 
          
    end
    disp('nv');
    disp(nv);
  //  if (nv)>1 then abort; nv=1; 
   // end
   // if (nv)<0 then abort; nv=0;
   // end
    disp(8); 
    nl=1-nv;
    x=0;
    y=0;
    for i=1:n
        x(i)=z(i)/(nl+nv*K(i));
        y(i)=x(i)*K(i);
    end
    disp('x'); disp(x);
    disp('y'); disp(y);
    disp('fun'); disp(fun);

endfunction
/*
O=[0.75630, 0.76974,  0.78017, 0.76921, 0.77698];
Zc=[0.33294, 0.31508,  0.30663, 0.31232, 0.31274];
Yi=[0.37447, 0.53248,  0.63875, 0.57554, 0.49550];
Tc=[190.56, 369.8,  408.14, 425.12, 305.32];
pc=[4599000, 4248000, 3648000,3769000,4872000];
w=[0.0115, 0.1523, 0.1770, 0.2002, 0.935];
z=[0.2,0.3,0.1,0.2,0.2];
*/
/*
O=[0.76974, 0.77698];
Zc=[0.31508,  0.31274];
Yi=[0.53248,  0.49550];
Tc=[369.8,  305.32];
pc=[4248000, 4872000];
w=[0.1523, 0.0935];
z=[0.5, 0.5];*/
/*
O=[ 0.78017];
Zc=[0.30663];
Yi=[ 0.63875];
Tc=[ 408.14];
pc=[ 3648000];
w=[0.1770];
z=[1];*/


O=[0.76974, 0.78017];
Zc=[0.31508,  0.30663];
Yi=[0.53248,  0.63875];
Tc=[369.8,  408.14];
pc=[4248000, 3648000];
w=[0.1523, 0.1770];
z=[0.8, 0.2];

/*
O=[0.77698, 0.76974, 0.76921];
Zc=[0.31274, 0.31508,  0.31232];
Yi=[0.49550, 0.53248,  0.44278];
Tc=[305.32, 369.8,  425];
pc=[4872000, 4248000, 4240000];
w=[0.0995, 0.1523, 0.2002];
z=[0.3, 0.4, 0.3];*/

a=[0,0,0,0,0];
b=[0,0,0,0,0];
c=[0,0,0,0,0];
d=[0,0,0,0,0];
/*p=680000;
T=260;*/
p=2000000;
T=350;

K=[0,0,0,0,0];
R=8.31;
n=2;
pt=0;
//n=2;
//O=[0.6974, 0. 75001];
//pc=[4248000, 2740000];
//Tc=[369.83, 540.2];
//w=[0.1523, 0.3495];
//Yi=[1.050+0.105*w(1)+0.482*w(1)^2,1.050+0.105*w(2)+0.482*w(2)^2];
//Zc=[0.31508, 0.3357-0.0294*w(2)];
//z=[0.7, 0.3];

for i=1:n
        K(i)=pc(i)*exp(5.373*(1+w(i))*(1-Tc(i)/T))/p;
        Kp(i)=pc(i)*exp(5.373*(1+w(i))*(1-Tc(i)/T))/p;
end
disp(K);

while 1==1 
disp(5);
for i=1:n
    al=O(i)^3;
    bet=Zc(i)+O(i)-1;
    sig=-Zc(i)+O(i)*(0.5+(O(i)-0.75)^0.5);
    del=-Zc(i)+O(i)*(0.5-(O(i)-0.75)^0.5);    
    ac=al*R^2*Tc(i)^2/pc(i);
    //a(i)=ac*Yi(i);
    a(i)=ac*(1+Yi(i)*(1-sqrt(Tc(i)/T)))^2;
    b(i)=bet*R*Tc(i)/pc(i);
    c(i)=sig*R*Tc(i)/pc(i);
    d(i)=del*R*Tc(i)/pc(i);

end

disp(3);
[x,y]=ras(z,K,n);
disp(6);
am=0;
for i=1:n
 //   am=am+y(i)*sqrt(a(i));
    for j=1:n
        am=am+y(i)*y(j)*sqrt(a(i)*a(j));            
    end
end
//am=am*am;
bm=0;

for i=1:n
    for j=1:n
        bm=bm+y(i)*y(j)*0.5*(b(i)+b(j)); 
    end
end

cm=0;
for i=1:n

        cm=cm+y(i)*c(i);

end

dm=0;
for i=1:n

        dm=dm+y(i)*d(i);

end

Am=am*p/(R^2*T^2);
Bm=bm*p/(R*T);
Cm=cm*p/(R*T);
Dm=dm*p/(R*T);
Bi=b.*p/(R*T);
Ci=c.*p/(R*T);
Di=d.*p/(R*T);

eq=(poly([-(Bm*Cm*Dm+Cm*Dm+Am*Bm) Am-Bm*Cm+Cm*Dm-Bm*Dm-Dm-Cm Cm+Dm-Bm-1 1] ,'x','c'));
ro=(roots(eq));
j=1;
z_g=0;
disp('ro zg');
disp(ro);

for i=1:3
    if (abs(imag(ro(i)))<1e-6) & (real(ro(i))>z_g)  then z_g=real(ro(i)); end;
end
disp(z_g); 
for i=1:n
        Am=real(Am);
        Bm=real(Bm);
        Cm=real(Cm);
        Dm=real(Dm);
        am=real(am);
        bm=real(bm);
        cm=real(cm);
        dm=real(dm);
        E1(i)=log(y(i)*p);
        E2(i)=log(z_g-Bm);
        
       
        aij=0;
        for j=1:n aij=aij+y(j)*sqrt(a(i)*a(j)); end;
        E3(i)=-Am/(Cm-Dm);
        E4(i)=(2*aij/am-(c(i)-d(i))/(cm-dm)); 
        E41(i)=log((z_g+Cm)/(z_g+Dm));
        E5(i)=+Bi(i)/(z_g-Bm);
        E6(i)=Am/(Cm-Dm)*(Ci(i)/(z_g+Cm)-Di(i)/(z_g+Dm));
        Fg(i)=(log(y(i)*p)-log(z_g-Bm)-Am/(Cm-Dm)*(2*aij/am-(c(i)-d(i))/(cm-dm))*log((z_g+Cm)/(z_g+Dm))+Bi(i)/(z_g-Bm)-Am/(Cm-Dm)*(Ci(i)/(z_g+Cm)-Di(i)/(z_g+Dm)));
        Fgp(i)=log(y(i)*p)-log(z_g-Bm)+b(i)/bm*(z_g-1)-Am/(2*sqrt(2)*Bm)*(2*aij/am-b(i)/bm)*log((z_g+(1+sqrt(2)*Bm)/(z_g+(1-sqrt(2)*Bm))));           
end
//break;    
/////////////////////////////////////////////////


am=0;
for i=1:n
 //   am=am+x(i)*sqrt(a(i));
    for j=1:n
        am=am+x(i)*x(j)*sqrt(a(i)*a(j));            
    end
end
//am=am*am;
bm=0;

for i=1:n
    for j=1:n
        bm=bm+x(i)*x(j)*0.5*(b(i)+b(j));            
    end /* bm=bm+x(i)*b(i);
        disp('x');disp(x(i));
        disp('b');disp(b(i));*/
end

cm=0;
for i=1:n

        cm=cm+x(i)*c(i);

end

dm=0;
for i=1:n

        dm=dm+x(i)*d(i);

end

Am=am*p/(R^2*T^2);
Bm=bm*p/(R*T);
Cm=cm*p/(R*T);
Dm=dm*p/(R*T);
Bi=b.*p/(R*T);
Ci=c.*p/(R*T);
Di=d.*p/(R*T);

eq=(poly([-(Bm*Cm*Dm+Cm*Dm+Am*Bm) Am-Bm*Cm+Cm*Dm-Bm*Dm-Dm-Cm Cm+Dm-Bm-1 1] ,'x','c'));
ro=(roots(eq));
disp('ro zl');
disp(ro);
j=1;
z_l=10000;

for i=1:3
    if (abs(imag(ro(i)))<1e-6) & (real(ro(i))<z_l) & (real(ro(i))>0) then z_l=real(ro(i)); end;
end
disp(z_l);
for i=1:n
        Am=(Am);
        Bm=(Bm);
        Cm=(Cm);
        Dm=(Dm);
        am=(am);
        bm=(bm);
        cm=(cm);
        dm=(dm);
        E1(i)=log(x(i)*p);
        E2(i)=log(z_l-Bm);
        
       
        aij=0;
        for j=1:n aij=aij+x(j)*sqrt(a(i)*a(j)); end;
        E3(i)=-Am/(Cm-Dm);
        E4(i)=(2*aij/am-(c(i)-d(i))/(cm-dm)); 
        E41(i)=log((z_l+Cm)/(z_l+Dm));
        E5(i)=+Bi(i)/(z_l-Bm);
        E6(i)=Am/(Cm-Dm)*(Ci(i)/(z_l+Cm)-Di(i)/(z_l+Dm));
        Fl(i)=(log(x(i)*p)-log(z_l-Bm)-Am/(Cm-Dm)*(2*aij/am-(c(i)-d(i))/(cm-dm))*log((z_l+Cm)/(z_l+Dm))+Bi(i)/(z_l-Bm)-Am/(Cm-Dm)*(Ci(i)/(z_l+Cm)-Di(i)/(z_l+Dm)));
        Flp(i)=log(x(i)*p)-log(z_l-Bm)+b(i)/bm*(z_l-1)-Am/(2*sqrt(2)*Bm)*(2*aij/am-b(i)/bm)*log((z_l+(1+sqrt(2)*Bm)/(z_l+(1-sqrt(2)*Bm))));
end; 
 //break;  
 //Fl=real(Fl);
pl=pl+1;
PL(pl)=K(1);
disp('FL');
disp(Fl); //disp(Flp);
disp('FG');
disp(Fg); //disp(Fgp);
pt=pt+1;
for i=1:n 
K1(pt,i)=K(i);
Kp1(pt,i)=Kp(i);
//K2(pt)=K(2);
end;

for i=1:n 
    K(i)=K(i)*Fl(i)/Fg(i);//exp(Fl(i))/exp(Fg(i));
    Kp(i)=Kp(i)*Fl(i)/Fg(i);
end
disp(K);


//break;
if pt>1 then ko=0 else ko=1 end;

if pt>1 then
    
    disp(10);
    /*disp(abs(K1(pt)/K1(pt-1)-1));
    disp(abs(K2(pt)/K2(pt-1)-1));
    disp(11);
    if (abs(K1(pt)/K1(pt-1)-1)<=0.1) & (abs(K2(pt)/K2(pt-1)-1)<=0.1) then break; end; */
//    for i=1:n
//    if (abs(K1(pt,i)/K1(pt-1,i)-1)>0.1) then ko=1; end;
//    end

for i=1:n
    disp(10);
    disp(abs(Fl(i)/Fg(i)-1));
    if abs(Fg(i)/Fl(i)-1)>0.001 then ko=1; end;
    disp(11);
end

end
if ko==0 then break; end;
disp(4);
 
//if  abs(K1(pt))>1000 then break; end;
end
