function f = momicDerivative6(obj, t, logintMu, shearRate)
% Derivatives of the 7 integer moments to be used in method of moments with
% interpolative closure (MOMIC).

intMu=exp(logintMu); 
f=zeros(6,1);

%Model Parameters
a_p   = obj.cnst.a_p;    % primary particle radius
phi_p = obj.cnst.phi_p;  % particle volume fraction
k_b   = obj.cnst.k_b;    % Boltzmann constant
T     = obj.cnst.T;      % Temperature
mu_s  = obj.cnst.mu_s;   % Medium viscosity (Pa-s)

df    = obj.par.d_f;     % Fractal dimension     
alfa  = obj.par.alfa;    % Collision efficiency
W     = obj.par.W;       % Fuch's stability ratio
b_o   = obj.par.b_0;     % Breakage constant
m_p   = obj.par.m_p;     % Number of primary particles in a primary cluster

fra_moment = @(j, c) obj.fra_moment(j, c);
theta = @(x) obj.theta(x);
shearcoe=phi_p*alfa*shearRate*obj.cutOff(logintMu)/(2*pi);

browncoe=phi_p*k_b*T*obj.cutOff(logintMu)...
    /(4*pi*(a_p^3)*mu_s*obj.etaTrimodal(logintMu)*W);

breakcoe=(shearRate^2)*b_o;

c=obj.MOMIC(logintMu);
 % 0th moment 
    shearAggregation0= -2*shearcoe*...
        (intMu(1)*fra_moment(3/df,c)+...
        3*fra_moment(2/df,c)*fra_moment(1/df,c));
    
    brownAggregation0= -2*browncoe*...
        (intMu(1)^2+fra_moment(1/df,c)*fra_moment(-1/df,c));
    
    breakAge0= breakcoe*(theta(0)-1)*...
        (fra_moment(1/df,c)-m_p^(1/df)*intMu(1));
    
    f(1)= (shearAggregation0 + brownAggregation0 + breakAge0)/intMu(1);
 % 2nd moment
    shearAggregation2= shearcoe*...
        (4*1*fra_moment(1+3/df,c)+...
        12*fra_moment(1+2/df,c)*fra_moment(1+1/df,c));
    %checked
    brownAggregation2=browncoe*4*(1^2+...
        fra_moment(1+1/df,c)*fra_moment(1-1/df,c));
    %checked  
    breakAge2= breakcoe*(theta(2)-1)*...
        (fra_moment(2+1/df,c)-m_p^(1/df)*intMu(2));
    
    f(2)= (shearAggregation2 +brownAggregation2 + breakAge2)/intMu(2);
           
% 3rd moment
    shearAggregation3= shearcoe*2*3*...
        (1*fra_moment(2+3/df,c)+...
        intMu(2)*fra_moment(1+3/df,c)+...
        3*fra_moment(2+2/df,c)*fra_moment(1+1/df,c)+...
        3*fra_moment(1+2/df,c)*fra_moment(2+1/df,c));
    %checked
    brownAggregation3= browncoe*2*3*...
        (2*1*intMu(2)+...
        fra_moment(2+1/df,c)*fra_moment(1-1/df,c)+...
        fra_moment(2-1/df,c)*fra_moment(1+1/df,c));
    %checked
    breakAge3= breakcoe*(theta(3)-1)*...
        (fra_moment(3+1/df,c)-m_p^(1/df)*intMu(3));
    
    f(3)= (shearAggregation3 + brownAggregation3 + breakAge3)/intMu(3);

%4th moment
    shearAggregation4= shearcoe*...
        (2*4*(intMu(3)*fra_moment(1+3/df,c)+...
        1*fra_moment(3+3/df,c)+...
        3*fra_moment(1+2/df,c)*fra_moment(3+1/df,c)+...
        3*fra_moment(1+1/df,c)*fra_moment(3+2/df,c))+...
        2*6*(intMu(2)*fra_moment(2+3/df,c)+...
        3*fra_moment(2+2/df,c)*fra_moment(2+1/df,c)));
    %checked
    brownAggregation4= browncoe*...
        (6*(2*intMu(2)^2+...
        2*fra_moment(2+1/df,c)*fra_moment(2-1/df,c))+...
        2*4*(2*intMu(3)*1+...
        fra_moment(3+1/df,c)*fra_moment(1-1/df,c)+...
        fra_moment(3-1/df,c)*fra_moment(1+1/df,c)));
    %checked
    breakAge4= breakcoe*(theta(4)-1)*...
        (fra_moment(4+1/df,c)-m_p^(1/df)*intMu(4));
    
    f(4)= (shearAggregation4 + brownAggregation4 + breakAge4)/intMu(4);
    
% 5th moment
    shearAggregation5= shearcoe*2*...
        (5*(1*fra_moment(4+3/df,c)+...
        3*fra_moment(4+2/df,c)*fra_moment(1+1/df,c)+...
        3*fra_moment(4+1/df,c)*fra_moment(1+2/df,c)+...
        intMu(4)*fra_moment(1+3/df,c))+...
        10*(fra_moment(3+3/df,c)*intMu(2)+...
        3*fra_moment(3+2/df,c)*fra_moment(2+1/df,c)+...
        3*fra_moment(2+2/df,c)*fra_moment(3+1/df,c)+...
        fra_moment(2+3/df,c)*intMu(3)));
   %checked 
    brownAggregation5= browncoe*...
        (2*5*(2*1*intMu(4)+...
        fra_moment(4+1/df,c)*fra_moment(1-1/df,c)+...
        fra_moment(4-1/df,c)*fra_moment(1+1/df,c))+...
        2*10*(2*intMu(3)*intMu(2)+...
        fra_moment(3+1/df,c)*fra_moment(2-1/df,c)+...
        fra_moment(3-1/df,c)*fra_moment(2+1/df,c)));
    %checked
    breakAge5= breakcoe*(theta(5)-1)*...
        (fra_moment(5+1/df,c)-m_p^(1/df)*intMu(5));
    
    f(5)= (shearAggregation5 + brownAggregation5 + breakAge5)/intMu(5);
    
% 6th moment
    shearAggregation6= shearcoe*...
        (2*6*(fra_moment(5+3/df,c)*1+...
        3*fra_moment(5+2/df,c)*fra_moment(1+1/df,c)+...
        3*fra_moment(5+1/df,c)*fra_moment(1+2/df,c)+...
        intMu(5)*fra_moment(1+3/df,c))+...
        2*15*(fra_moment(4+3/df,c)*intMu(2)+...
        3*fra_moment(4+2/df,c)*fra_moment(2+1/df,c)+...
        3*fra_moment(4+1/df,c)*fra_moment(2+2/df,c)+...
        intMu(4)*fra_moment(2+3/df,c))+...
        2*20*(fra_moment(3+3/df,c)*intMu(3)+...
        3*fra_moment(3+2/df,c)*fra_moment(3+1/df,c)));
    %checked
    brownAggregation6= 2*browncoe*...
        (2*6*(1*intMu(5)+...
        fra_moment(5+1/df,c)*fra_moment(1-1/df,c)+...
        fra_moment(5-1/df,c)*fra_moment(1+1/df,c))+...
        2*15*(intMu(4)*intMu(2)+...
        fra_moment(4+1/df,c)*fra_moment(2-1/df,c)+...
        fra_moment(4-1/df,c)*fra_moment(2+1/df,c))+...
        2*20*(intMu(3)^2+...
        fra_moment(3+1/df,c)*fra_moment(3-1/df,c)));
    %checked
    breakAge6= breakcoe*(theta(6)-1)*...
        (fra_moment(6+1/df,c)-m_p^(1/df)*intMu(6));
    
    f(6)= (shearAggregation6 + brownAggregation6 + breakAge6)/intMu(6);
    
end