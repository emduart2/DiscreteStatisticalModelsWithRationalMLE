--- In this script we define the functions that take care of the computations
--- in "Discrete Statistical Models with Rational MLE."
---

--------------------------------------------------
--------------------------------------------------
--*       hUniform     *--------------------------
--------------------------------------------------
-- This functions takes as input a matrix B and has 
-- output the factors of the psi coordinates of the
-- Horn uniformization corresponding to B.
-- each element in the list is a list of the factors
-- of the coordinate functions. These are stored
-- as Expressions of type Product
hUniform=(B)->(
    numberOfAtoms= rank source B;
    uVector=apply(numberOfAtoms, i-> u_i);
    Rr= QQ[uVector];
    uVector=apply(numberOfAtoms, i-> u_i);
    uVector=transpose matrix{uVector};
    uB=B*uVector;
    S=frac Rr;
    uB=sub(uB,S);
    j=0;
    unif={};
    num={};
    den={};
    while j< rank source B do(
    factors=for i from 0 when i < rank target uB list (if (B_(i,j) == 0) then continue 1; Power{uB_(i,0),(B_(i,j))});
    factors=delete(1,factors);
    unif=unif|{factors};
    j=j+1;
    );
    return unif
)

--------------------------------------------------
--------------------------------------------------
--*       toFraction     *------------------------
--------------------------------------------------
-- This functions takes a list whose elements are
-- expressions of class product and decides
-- whether the element is a numerator or a denominator
-- it returns a fraction which is how we would write
-- a given psi coordinate
toFraction=(factors)->(
    i=0;
    num={};
    den={};
    while i< #factors do(
     if factors#i#1>0 then(
	 num=num|{factors#i};
	 --print num;
	 );
     if factors#i#1<0 then(
	 den=den|{Power{factors#i#0,-factors#i#1}};
	 --print den;
	 );
     i=i+1;
     );
     quot=Divide{Product num,Product den};
     return quot
    )

--------------------------------------------------
--------------------------------------------------
--*       psiCoord     *--------------------------
--------------------------------------------------
-- This function takes the matrix B and returns
-- the psi coordinates of the function already written
-- as nice fractions.
psiCoord=(B)->(
    uniB=hUniform B;
    psi= apply(uniB,i-> toFraction i);
    return psi
    )
--------------------------------------------------
--------------------------------------------------

--------------------------------------------------
--------------------------------------------------
--*       checkFriendly     *---------------------
--------------------------------------------------
-- This function takes the matrix B and returns
-- either the vector lambda that makes a matrix friendly
-- or a satement saying that no such lambda exists.
--------------------------------------------------
--------------------------------------------------

checkFriendly=(B)->(
    psiFactors= apply(hUniform B,i-> product i);
    psiNumerators= apply(psiFactors,i-> numerator value i);
    psiDenominators=apply(psiFactors,i-> denominator value i);
    n=#psiFactors;
    lambdas=apply(n,i->l_i);
    LL=QQ[flatten entries uVector,lambdas];
    L=frac(LL)    ;
    lambda=apply(n,i->l_i);
    psiLambdaNum=apply(#lambda,i-> (sub(psiNumerators#i,L)*lambda#i));
    use L;
    psiLambdaSum=sum apply(#psiLambdaNum,i->psiLambdaNum#i/sub(psiDenominators#i,L));
    num = numerator psiLambdaSum;
    den = denominator psiLambdaSum;
    (degree den)#0 +1 ==    (degree num)#0; -- if these degrees don't match then there is no chance.
    basisU=super basis((degree den)#0,Rr);
    use LL;
    lambdaEquations = contract(sub(basisU,LL),num);
    lambdaConstants = contract(sub(basisU,LL),den);
    lambdaEquations_(0,0);
    rowsCoefficientMatrix = apply(flatten entries lambdaEquations,
	i->contract(sub(matrix{lambda},LL),i));
    CoefficientMatrix = fold(rowsCoefficientMatrix,(i,j)-> i||j);
    if rank(CoefficientMatrix|(transpose lambdaConstants))> rank CoefficientMatrix then return "No lambda exists"-- If true, then lambda doesn't exist
    else return lambdaSolution = solve(CoefficientMatrix,transpose lambdaConstants)
)
--------------------------------------------------
--------------------------------------------------
--*       checkPair     *---------------------
--------------------------------------------------

-- Check sum to one conditions on a matrix B
-- the code for checkFriendly takes as input a pair (B,lambda)
-- and checks if the given pair is a friendly pair.
--------------------------------------------------
--------------------------------------------------
checkPair=(B,Lambda)->(
    use S;
    psiFactors= apply(hUniform B,i-> product i);
    psiFactors= matrix{apply(psiFactors,i-> value i)};
    sumPsi = psiFactors *(transpose sub(Lambda,S));
    return sumPsi
    )
--------------------------------------------------
--------------------------------------------------

--------------------------------------------------
--------------------------------------------------
--*       implicitEquations     *---------------------
--------------------------------------------------

-- Uses an elimination algorithm for rational functions
-- to compute the implicit equations of a rational 
-- map.
--------------------------------------------------
--------------------------------------------------


implicitEquations =  (Bmat,lam)->(
     psiC := psiCoord Bmat; 
     n  := #psiC;
     uVector := apply(n, i-> u_i);
     Distr=QQ[y,z,p_0..p_(n-1),uVector];
     uVector=apply(n, i-> u_i);
     uVector=transpose matrix{uVector};
     --if (checkFriendly Bmat === "No lambda exists") then ( return  "No lambda exists")
     --else(
	 lambdaMultiplier = flatten entries sub(lam,Distr);
          f=apply(n,i->value( sub((denominator psiC_i),Distr)*p_i- lambdaMultiplier_i*z*sub(numerator psiC_i,Distr)) );
          lcmden = lcm apply(psiC, i-> sub(value denominator i,Distr));
          I=ideal( f|{1-lcmden*y});
          return eliminate({z}|flatten entries sub(uVector,Distr),I)--)
     )
     

--------------------------------------------------

--------------------------------------------------

--------------------------------------------------
--------------------------------------------------
--*       pertinentTerm     *--------------------------
--------------------------------------------------
-- This function takes a polynomial Da and a term in in Da
-- and checks whether the system of equations G_l has a solution
-- It returns Not pertinent if the system G_l does not have a solution.
--  It returns the solution if the system has one.
-- ALGORITHM: We make a list of exponents of monomials in Da and
-- delete from it the element in the list that corresponds to the
-- monomial termDa. To the each exponenvector in this list
-- we substract the exponent vector of termDa. This gives us a list
-- of the exponents of the monomials in g=1-\sum c_lambda*monomials
-- We call this lis exponentlist, this is our matrix B corresponding
-- to the pair(Da,termDa). Next we take the matrix B modulo 2 and tranpose
-- it to get the coefficients of the system.
-- Now we extract the vector of coefficients of Da using newExponents and
-- cLambda. Afterwards we compute cSign to keep track using 0,1 por positive and
-- negative coefficients respectively.
-- If the system does not have a solution we return "Not pertinent"
-- If the system does have a solution we return the the set of solutions.
-- IT could happen that there are feasible solutions in the system of equations
-- that we set up but when we evaluate at a positive number the image is negative, therefore we
-- need to incorporate the positivitycheck
pertinentTerm=(Da,termDa)->(
    currentMonomial := termDa;
    currentC := coefficient(toHMonomial(flatten exponents currentMonomial),Da);
    exponentsDa := exponents Da;
    zeros := for i from 1 to  (numgens DualCoord) list 0;
    i:=-1;
    Bcolumns :={};
    cLambda :={};
    while  (i<#exponentsDa-1) do(
       i=i+1;
       newExponent = exponentsDa#i-(flatten exponents currentMonomial);
       if (newExponent=== zeros) then continue;
       newC     = -coefficient(toHMonomial exponentsDa#i,Da)/currentC;
       Bcolumns = Bcolumns|{newExponent};
       cLambda  = cLambda|{newC};
       );
    B :=  transpose matrix Bcolumns;
    Bmod2 := sub(B, ZZ/2);
    Bmod2t  := transpose Bmod2;
    cSign:=transpose matrix{for i from 0 to #cLambda-1 list( if  (cLambda#i>=0)  then continue 0;1)}**(ZZ/2);
    if rank(Bmod2t)== rank(Bmod2t|cSign) then return   (B,cLambda,currentMonomial,positivityCheck(B,cLambda))--((gens ker Bmod2t),solve(Bmod2t,cSign))
    else return "Not pertinent"
   )
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------


--------------------------------------------------
--------------------------------------------------
--*       pertinentDiscriminant     *-------------
--------------------------------------------------
--- INPUT: a discriminant Da in the ring DualCoord
--- OUTPUT: a list where each entry in the list 
--- to the result of testing whether a given term
--- in Da is pertinent.
pertinentDiscriminant=(Da)->(
    termsDa:= terms Da;
    return for i from 0 to #termsDa-1 list pertinentTerm(Da,termsDa_i)
    )
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------

--------------------------------------------------
--------------------------------------------------
--*       pertinentB               *--------------
--------------------------------------------------
-- This function takes a matrix B, figures out if it is
-- friendly or not and then decides if the pair is
-- pertinent or not based on finding solutions to
-- the system of equations
--------------------------------------------------
pertinentB=(B,lam)->(
    kk=ZZ/2;
    Bmod2t= transpose(sub(B,kk));
    --cLambda=flatten entries checkFriendly B;
    cLambda=flatten entries lam;
    cSign=transpose( matrix{for i from 0 to #cLambda-1 list( if  (sub(cLambda#i,QQ))>=0  then continue 0;1)})**kk;
    if rank(Bmod2t)== rank(Bmod2t|cSign) then return   ((gens ker Bmod2t),solve(Bmod2t,cSign))
    else return "Not pertinent"
    )

--------------------------------------------------
--------------------------------------------------
--------------------------------------------------



-----------------------------------------------------
--*    toMonomial *----------------------------------
-----------------------------------------------------
-- This functions takes a vector of entries and outputs
-- a monomial correspondic to x^vec.
-- The variables are hard wired to use with caution.
-----------------------------------------------------
toMonomial=(vec)->(
    variables:= flatten entries vars DualCoord;
    factors:= apply(#vec,i->(variables_i)^(vec#i));
    return fold(factors,(i,j)->i*j)
    )
-----------------------------------------------------
-----------------------------------------------------
-----------------------------------------------------

--- this function takes an A-discriminat Da and a term in it termDa
--- INPUT: 
--- Da= discriminant in the polynomial ring DualCoord
--- termDa= a term of Da
--- OUTPUT:
--- A friendly pair (H,lambda)
--- 
--- and outputs the Bhorn matrix obtained after factoring the termDa
-- from Da, i.e. termDa*(1-\sum c_lambda_j x^b_j) 
-- The function returns

hornTermDa=(Da,termDa)->(
   currentMonomial := termDa;
   currentC := coefficient(toHMonomial(flatten exponents currentMonomial),Da);
   exponentsDa := exponents Da;
   zeros := for i from 1 to  (numgens DualCoord) list 0;
   i:=-1;
   Bcolumns :={};
   cLambda :={};
   while  (i<#exponentsDa-1) do(
       i=i+1;
       newExponent = exponentsDa#i-(flatten exponents currentMonomial);
       if (newExponent=== zeros) then continue;
       newC     = -coefficient(toHMonomial exponentsDa#i,Da)/currentC;
       Bcolumns = Bcolumns|{newExponent};
       cLambda  = cLambda|{newC};
       );
   B:= transpose matrix Bcolumns;
   c:= matrix{cLambda};
   return {B,c}
   )
-----------------------------------------------------
-----------------------------------------------------
-----------------------------------------------------

toHMonomial=(vec)->(
    variables := flatten entries vars DualCoord;
    factors := apply(#vec,i->(variables_i)^(vec#i));
    return fold(factors,(i,j)->i*j)
    )

-----------------------------------------------------
-----------------------------------------------------
--- Functions to compute affine A-discriminants
--- INPUT: a matrix A
--- OUTPUT: a discriminant Da.

aDiscriminantA=(A)->(
   r  :=rank target A;
   k  :=rank source A;
   Tt :=QQ[x_1..x_k,t_1..t_r];
   varsTt :=matrix{apply(r,i->t_(i+1))};
   f :=sum apply(k,i->(x_(i+1))*toAdiscrMonomial(r,flatten entries A_i));
   I := ideal(f,diff(varsTt,f));
   discrIdeal := eliminate(flatten entries varsTt,I);
   return  (gens discrIdeal)_(0,0)
   )

toAdiscrMonomial=(r,vec)->(
    tvars := apply(r,i->t_(i+1));
    factors := apply(#vec,i->(tvars_i)^(vec#i));
    return fold(factors,(i,j)->i*j)
    )
-----------------------------------------------------
--- This basic function takes a two vectors of integers 
--- and outputs a list of element i in vec to the power of element
--- i in the list exponentsVec
exponentiateVector = (vec,exponentsVec)->(
    factors:=apply(#vec,i->(sub(vec#i,RR))^(exponentsVec#i));
    return factors
    )
--- This function checks wether a given pair (H,lambda)
--- satisfies the positivity conditions for the image of a point 
--- to be positive.
positivityCheck = (H,lambda)->(
    rowSumH  := apply(entries H, i->sum i);
    exponentVectors :=entries transpose H;
    listCoordinates :=apply(#exponentVectors-1,i->sub(lambda#i,RR)*product(exponentiateVector(rowSumH,exponentVectors#i)));
    return all(listCoordinates,j->(j>0));
    )
    --return  all(apply(#exponentVectors,i->sub(-lambda#i,RR)*product(exponentiateVector(rowSumH,exponentVectors#i))),i->(i>0));
    
-----------------------------------------------------


end--

restart
load "~/Documents/gitHubRepos/DiscreteStatisticalModelsWithRationalMLE/functionsMLE.m2"
