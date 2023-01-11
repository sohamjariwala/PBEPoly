function qoptim = ParTemp(objfunc,q1,EBHot,constraint1,type1,constraint2,type2)
% Parallel temper algorithm to determine minimum of a user specified
% objective function, objfunc.  q1 is a vector containing the initial guess
% for the parameters involved in the objective function.  The user should
% also specify a hot energy level.  If the user fails to specify a level, 
% the hot energy level is assumed to be unitary.  This should be the same 
% order of magnitude as the objective function value.  The user may also 
% specify constraints on the parameters where the constraint value is a 
% vector of either maximum or minimum values for each paramters, and the 
% type is a string consisting of either 'max' or 'min'.

format long

% Check constraints on variables
type=0;
if nargin>5         %Two constraints provided by user
    if strcmpi(type1,'MAX')&&strcmpi(type2,'MIN')
        maxvals=constraint1;
        minvals=constraint2;
    else
        maxvals=constraint2;
        minvals=constraint1;
    end
    q1=q1-minvals;                  %Min constraint on variables
    q1=q1./(maxvals-minvals-q1);    %Max constraint on variables
    type=1;
elseif nargin>3     %One constraint provided by user
    if strcmpi(type1,'MAX')         %Max constraint provided
        maxvals=constraint1;
        q1=q1./(maxvals-q1);
        type=2;
    else                            %Min constraint provided
        minvals=constraint1;
        q1=q1-minvals;
        type=3;
    end
elseif nargin<3
    EBHot=1;
end

% Initialize parameters
Nrun=15;                %Number of parallel runs, change as needed
EBCold=EBHot*10^(-10);  %Cold energy level, change multiplier as needed
p=3;                    %MC standard deviation exponent, change as needed
Nmin=5;                 %Multiple of maximum NEx for termination criteria, change as needed
N=length(q1);           %Number of parameters
Rkcrit=0.5;             %Critical autocorrelation function value
A=ones(N,1);

% Calculate EB for each run
EB=zeros(Nrun,1);
EB(1)=EBHot;
for j=1:Nrun-1
    EB(j+1)=EB(j)*(EBCold/EBHot)^(1/(Nrun-1));
end

% Initial run
Fobj=callobj(q1);

% Initialization of loop variables for determination of NExHot and NExCold
theta=zeros(N,1);                   %Log parameter vector
sigma=zeros(N,1);                   %MC standard deviation
Nsteps=1000;                        %Autocorrelation number of steps
FobjHC=zeros(2,2*Nsteps);           %Objective function array
FobjHC(:,1)=Fobj;
EBHC=[EBHot,EBCold];                %Energy values vector

% Initialize parameter arrays
qHC=zeros(2,N,Nsteps*2);
qHC(1,:,1)=q1;
qHC(2,:,1)=q1;

fprintf('Preparing for Iterations...\n\n')
% Iterations without exchanges for Hot and Cold energy levels
for k=2:2*Nsteps
    for j=1:2
        
        % MC Step
        for n=1:N
            theta(n)=log(qHC(j,n,k-1));
            sigma(n)=(EBHC(j)/EBHot)^(1/p)*A(n);
            zeta=normrnd(0,1);
            theta(n)=theta(n)+sigma(n)*zeta;
            qHC(j,n,k)=exp(theta(n));
        end
        
        % Calculate new values for objective function
        if isrow(q1)
            FobjHC(j,k)=callobj(qHC(j,:,k));
        else
            FobjHC(j,k)=callobj(qHC(j,:,k)');
        end
        
        % Check acceptance of new values
        if FobjHC(j,k)<=FobjHC(j,k-1)
            Paccept=1;
        else
            Paccept=exp(-(FobjHC(j,k)-FobjHC(j,k-1))/EBHC(j));
        end
        if rand>Paccept
            FobjHC(j,k)=FobjHC(j,k-1);
            qHC(j,:,k)=qHC(j,:,k-1);
        end
    end
end

%Calculate autocorrelation function
acfHot=autocorr(qHC(1,1,:),Nsteps);
acfCold=autocorr(qHC(2,1,:),Nsteps);

%Determine NExHot and NExCold
NExHot=1;
while acfHot(NExHot)>Rkcrit && NExHot<Nsteps
    NExHot=NExHot+1;
end
NExCold=1;
while acfCold(NExCold)>Rkcrit && NExCold<Nsteps
    NExCold=NExCold+1;
end

% Calculate NEx for each run
NEx=zeros(Nrun,1);
NEx(1)=NExHot;
for j=1:Nrun-1
    NEx(j+1)=NEx(j)*(NExCold/NExHot)^(1/(Nrun-1));
end

% Re-Initialization of loop variables
FobjNew=zeros(Nrun,1);              %New values for objective function
qNew=zeros(Nrun,N);                 %New parameter values
theta=zeros(N,1);                   %Log parameter vector
sigma=zeros(N,1);                   %MC standard deviation
k=0;                                %Loop variable
% kred=0;                             %Optional reducing counter to enable higher precision
FobjBest=zeros(Nmin*NExCold+1,1);   %Best objective function for each iteration
exc=ones(Nrun,1);                   %Exchange signal vector
FobjOld=ones(Nrun,1)*Fobj;          %Starting objective function values
qOld=zeros(Nrun,N);
for j=1:Nrun
    qOld(j,:)=q1;                   %Starting parameter values
end

fprintf('Preparation Complete\n\n')
display(0,Fobj,q1)

% Iterative step with exchanges
while k<=Nmin*NExCold||abs(FobjBest(1)-FobjBest(Nmin*NExCold+1))>EBCold
    k=k+1;
    
    % Enable to increase precision accurate answers
%     if k>=Nmin*NExCold&&k-kred>Nmin*NExCold/4&&abs(FobjBest(round(Nmin*NExCold/4))-FobjBest(Nmin*NExCold+1))<=EBCold
%         kred=k;
%         A=A*0.001;
%     end

    FobjBestrun=Inf;
    for j=1:Nrun
        
        % MC Step
        for n=1:N
            theta(n)=log(qOld(j,n));
            sigma(n)=(EB(j)/EBHot)^(1/p)*A(n);
            zeta=normrnd(0,1);
            theta(n)=theta(n)+sigma(n)*zeta;
            qNew(j,n)=exp(theta(n));
        end
        
        % Calculate new values for objective function
        if isrow(q1)
            FobjNew(j)=callobj(qNew(j,:));
        else
            FobjNew(j)=callobj(qNew(j,:)');
        end
        
        % Check acceptance of new values
        if FobjNew(j)<=FobjOld(j)
            Paccept=1;
        else
            Paccept=exp(-(FobjNew(j)-FobjOld(j))/EB(j));
        end
        if rand<=Paccept
            FobjOld(j)=FobjNew(j);
            qOld(j,:)=qNew(j,:);
        else
            FobjNew(j)=FobjOld(j);
            qNew(j,:)=qOld(j,:);
        end
        
        % Check exchanges
        if k>=NEx(j)*exc(j)&&j>1
            exc(j)=exc(j)+1;
            if FobjOld(j)>FobjOld(j-1)
                Paccept=1;
            else
                Paccept=exp((1/EB(j)-1/EB(j-1))*(FobjOld(j)-FobjOld(1)));
            end
            
            % Exchange occurs
            if rand<=Paccept
                FobjOld(j)=FobjOld(j-1);
                qOld(j,:)=qOld(j-1,:);
                FobjOld(j-1)=FobjNew(j);
                qOld(j-1,:)=qNew(j,:);
            end
        end
    end
    
    % Update best objective function value for the current iteration
    if FobjNew(j)<FobjBestrun
        FobjBestrun=FobjNew(j);
        index=j;
    end
    
    % Shift and update FobjBest
    if k<=Nmin*NExCold+1
        FobjBest(k)=FobjBestrun;
    else
        FobjBest(1:Nmin*NExCold)=FobjBest(2:Nmin*NExCold+1);
        FobjBest(Nmin*NExCold+1)=FobjBestrun;
    end
    if ~rem(k,1)
        display(k,FobjBestrun,qNew(index,:))
    end
end
if isrow(q1)
    qoptim=qNew(index,:);
else
    qoptim=qNew(index,:)';
end

% Undo correction to parameters for min/max values
if type==1
    qoptim=(qoptim.*maxvals+minvals)./(1+qoptim);
elseif type==2
    qoptim=qoptim.*maxvals./(1+qoptim);
elseif type==3
    qoptim=qoptim+minvals;
end

fprintf('Optimal Solution Found\n')
fprintf('Total Iterations = %d\n',k)
fprintf('Final Objective Value = %d\n',FobjBestrun)
if N>1
fprintf('Optimal Parameters = [')
for n=1:N-1
    fprintf('%d, ',qoptim(n))
end
fprintf('%d]\n\n',qoptim(N))
else
    fprintf('Optimal Parameter = %d\n\n',qoptim(1))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj=callobj(q)
        % Transorms variable then calls objective function.
        if type
            if type==1
                q=(q.*maxvals+minvals)./(1+q);
            elseif type==2
                q=q.*maxvals./(1+q);
            elseif type==3
                q=q+minvals;
            end
        end
        obj=objfunc(q);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function display(k,obj,q)
        % Transorms variable then displays iteration number, objective 
        % value, and current parameter values.
        if type
            if type==1
                q=(q.*maxvals+minvals)./(1+q);
            elseif type==2
                q=q.*maxvals./(1+q);
            elseif type==3
                q=q+minvals;
            end
        end
        fprintf('Iteration #%d\n',k)
        fprintf('Objective = %d\n',obj)
        if N>1
            fprintf('Current Parameters = [')
            for ndisp=1:N-1
                fprintf('%d, ',q(ndisp))
            end
            fprintf('%d]\n\n',q(N))
        else
            fprintf('Current Parameter = %d\n\n',q(1))
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end