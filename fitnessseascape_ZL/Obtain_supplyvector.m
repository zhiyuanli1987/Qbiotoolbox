function [FBv_givenIa,onelargesupplypoint,fval]=Obtain_supplyvector(metabfun,env)

fval=0;
intercell_dim=metabfun.intercell_dim;
largev=10;

%the steady state xs under given environment
withincellsteady=nan;
if intercell_dim>0
    [growthrate,withincellsteady,fval]=Obtain_Growth_givenenv(env,metabfun);
end

Intakes=metabfun.Intakes;
intakerates=Intakes(env,withincellsteady);


if intakerates(1)==0
    FBv_givenIa=@(Ia)env+[0;largev];
    onelargesupplypoint=env+[0;largev];
else
    
    FBv_givenIa=@(Ia)env+intakerates/intakerates(1)*(Ia-env(1));
    if intakerates(1)<0
        onelargesupplypoint=FBv_givenIa(-1*largev);
    else
        onelargesupplypoint=FBv_givenIa(largev);
    end
    
end

