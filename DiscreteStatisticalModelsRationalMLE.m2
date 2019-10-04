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
-- friendly or nor and then decides if the pair is
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
--
--


restart
-- CODE INSTRUCTIONS:
-- Run the following command to load all the functions. Change the path according to where you have placed
-- the aDiscriminants.m2 files
-- List of functions = {hUniform,toFraction,psiCoord,checkFriendly,checkPair,implicitEquations,pertinentTerm,
-- ,pertinentDiscriminant,toMonomial,hornTermDa,toHmonomials,pertinentB}
load "~/Documents/gitHubRepos/DiscreteStatisticalModelsWithRationalMLE/DiscreteStatisticalModelsRationalMLE.m2"
-- start using the functions to  try out examples using your favorite B matrices
-- EXAMPLE 1 -- Basic usage

A         =  matrix{{0,1,2,3}}
DualCoord = QQ[x_1,x_2,x_3,x_4]
Da        =  sub(aDiscriminantA A, DualCoord)
pertinentDiscriminant(sub(Da,DualCoord))
hornTermDa(sub(Da,DualCoord),sub(term,DualCoord))
term=-27*x_1^2*x_4^2
pertinentTerm(sub(Da,DualCoord),sub(term,DualCoord))
pertinentTerm(sub(x_1+x_2+x_3+x_4,DualCoord), sub(x_3,DualCoord))

DualCoord = QQ[x_1,x_2,x_3,x_4]
Dacubic= x_2^2*x_3^2-4*x_1*x_3^3-4*x_2^3*x_4+18*x_1*x_2*x_3*x_4-27*x_1^2*x_4^2
pertinentDiscriminant(Dacubic)
tDacubic= terms Dacubic
pair=hornTermDa(Dacubic,tDacubic#4)
pertinentB(pair#0,pair#1)
positivityCheck(pair#0,flatten entries pair#1)

----------------------
-- Example 26. 
----------------------
--For distinct positive integers alpha,beta,gamma with gcd(alpha,beta,gamma)=1, let
-- A = (0,alpha,beta,gamma).
-- First we generate all the A-matrices for which we will calculate the discriminant.
-- We call the list of all such matrices Amatrices.

--- These lines list all tuples of the form  {0,alpha,beta,gamma}. This produces a list L.
i=1
L={}
while i<16 do(
    j=i+1;
    while j<17 do(
	k=j+1;
	while k<18 do(
        L=L|{{0,i,j,k}};
	k=k+1;
	);
        j=j+1;
    );
    i=i+1;
    )
---
--- We create the list of all allowed tuples, called Amatrices.
Amatrices = for i from 0 to #L-1 list( if(gcd(L#i)!=1) then continue; L#i)
#Amatrices  -- There are 613 matrices in Amatrices.
--- 
--- In the next step we compute the discriminant of all the entries in Amatrices.
--- This is long. The user can skip this step and use the saved list trinomialDiscriminants
--- It is important to specify what is the ring DualCoord in which the polynomials
--- in the list trinomialDiscrims live. This is important to
--- call the function pertinentDiscriminant.
DualCoord=QQ[x_1..x_4]
trinomialDiscrims = time apply(Amatrices,i->sub(aDiscriminantA(matrix{i}),DualCoord));
#trinomialDiscrims
"trinomialDiscrims" << toString trinomialDiscrims << endl << close -- saved output to file
-- Need to make sure to load the correct DualCoord ring so that the code works.
DualCoord=QQ[x_1..x_4]
testTrinomialDiscrims= for j from 0 to #trinomialDiscrims-1 list (pertinentDiscriminant(sub(trinomialDiscrims#j,DualCoord)));
rawResults = delete("Not pertinent",flatten testTrinomialDiscrims);
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
for i from 0 to 10 list statisticalModels#i
----------------------
-- Example 27. 
----------------------
i=1
V={}
while i< 21 do(
    j=i+1;
    while j< 22 do(
	V=V|{{0,i,j}};
	j=j+1;
	);
    i=i+1;
    )
V
Vtotal=for i from 0 to #V-1 list( if(gcd(V#i)!=1) then continue; V#i);
i=0
j=0
Vpairs={}
while i<#Vtotal-1 do(
    while j<#Vtotal-1 do(
	Vpairs=Vpairs|{{Vtotal#i|Vtotal#j}};
	j=j+1;
	);
    i=i+1;
    )
Vpairs#0
#Vpairs
DualCoord = QQ[x_1..x_6]
resultantDiscrims =  time for k from 0  to #Vpairs-1 list(sub(aDiscriminantA((matrix Vpairs#k)||matrix{{0,0,0,1,1,1}}),DualCoord));
#resultantDiscrims

DualCoord = QQ[x_1..x_6] -- Make sure the ring is the correct one.
Da=sub(resultantDiscrims#0,DualCoord)
pertinentDiscriminant(Da)
testResultantDiscrims = for j from 0 to #resultantDiscrims-1 list (pertinentDiscriminant(sub(resultantDiscrims#j,DualCoord)));
rawResults = delete("Not pertinent",flatten testResultantDiscrims);
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
"trinomialDiscrims" << toString trinomialDiscrims << endl << close -- saved output to file
----------------------
-- Example 28. 
----------------------

DualCoord=QQ[x_1..x_4]
DA = x_1+x_2+x_3+x_4
--- Code for postive binomial multiples
a=1
deg=9
binomialPlusMultiples={};
while a < deg do(
bzz= flatten entries super basis(a,DualCoord);
binomials = subsets(bzz,2);
currentBinomialMultiple = for j from 0 to #binomials-1 list( if(gcd(binomials#j)!=1) then continue ; DA*(sum binomials#j));
print(#currentBinomialMultiple);
binomialPlusMultiples=binomialPlusMultiples|currentBinomialMultiple;
a=a+1;
)
#binomialPlusMultiples
binomialPlusMultiples#0
testBinomialPlusMultiples = for j from 0 to #binomialPlusMultiples -1 list pertinentDiscriminant(binomialPlusMultiples#j);
#(flatten testBinomialPlusMultiples)
rawResults = delete("Not pertinent",flatten testBinomialPlusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
statisticalModels
"stModelsBinomialPlusMultiples" << toString statisticalModels << endl << close -- saved output to file
--- Code for negative binomial multiples
a=1
deg=9
binomialMinusMultiples={};
while a < deg do(
bzz= flatten entries super basis(a,DualCoord);
binomials = subsets(bzz,2);
currentBinomialMultiple = for j from 0 to #binomials-1 list( if(gcd(binomials#j)!=1) then continue ; DA*(sum {(binomials#j)#0,-1*(binomials#j)#1}));
print(#currentBinomialMultiple);
binomialMinusMultiples=binomialMinusMultiples|currentBinomialMultiple;
a=a+1;
)
#binomialMinusMultiples
testBinomialMinusMultiples = for j from 0 to #binomialPlusMultiples -1 list pertinentDiscriminant(binomialMinusMultiples#j);
#(flatten testBinomialMinusMultiples)
rawResults = delete("Not pertinent",flatten testBinomialMinusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
statisticalModels
---------------------------------------------
--- Code for positive trinomial multiples
---------------------------------------------
a=1
deg=4
trinomialPlusMultiples={};
while a < deg do(
bzz= flatten entries super basis(a,DualCoord);
trinomials = subsets(bzz,3);
currentBinomialMultiple = for j from 0 to #trinomials-1 list( if(gcd(trinomials#j)!=1) then continue ; DA * (sum trinomials#j));
print(#currentBinomialMultiple);
trinomialPlusMultiples=trinomialPlusMultiples|currentBinomialMultiple;
a=a+1;
)
#trinomialPlusMultiples
testTrinomialPlusMultiples = for j from 0 to #trinomialPlusMultiples -1 list pertinentDiscriminant(trinomialPlusMultiples#j);
#(flatten testTrinomialPlusMultiples)
rawResults = delete("Not pertinent",flatten testTrinomialPlusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
statisticalModels
---------------------------------------------
--- Code for negative trinomial multiples
---------------------------------------------
a=1
deg=4
trinomialMinusMultiples={};
while a < deg do(
bzz= flatten entries super basis(a,DualCoord);
trinomials = subsets(bzz,3);
currentBinomialMultiple = for j from 0 to #trinomials-1 list( if(gcd(trinomials#j)!=1) then continue ; DA * (sum {(trinomials#j)#0,(trinomials#j)#1,-1*(trinomials#j)#2}));
print(#currentBinomialMultiple);
trinomialMinusMultiples=trinomialMinusMultiples|currentBinomialMultiple;
a=a+1;
)
#trinomialMinusMultiples
testTrinomialMinusMultiples = for j from 0 to #trinomialMinusMultiples -1 list pertinentDiscriminant(trinomialMinusMultiples#j);
#(flatten testTrinomialMinusMultiples)
rawResults = delete("Not pertinent",flatten testTrinomialMinusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
statisticalModels
--------


checkIrreducible = (form)->(
    Factors := factor form;
    if (#Factors) != 1 then return false; 
    return true;
    )
checkIrreducible(x_1^2-x_2^2)

#binomialPlusMultiples




binomialMultiples
\subsets(oo,2)
for j from 0 to #oo-1 list( if(gcd(oo#j)!=1) then continue; oo#j)
statisticalModels#0
last rawResults#7
for j from 0 to 10 list (if j==2 then continue; if j== 3 then continue;j)
instance(o56#1,String)

delete(false,noPertinent)

first (1,2)

aRawResults=for i from 0 to #rawResults-1 list (if (rawResults#i=="Not pertinent") then continue; i)

flatten o133
#o136
pertinentDiscriminant(sub(trinomialDiscrims#3,DualCoord))

B=matrix{{1,0,0},{0,1,0},{0,0,1},{-1,-1,-1}}
pertinentB(B,matrix{{1,1,1}})
B=matrix{{2,1,0},{0,1,1},{-2,-2,-1}}; lam={1,1,-1}
pertinentB(B,matrix{lam})
DualCoord=QQ[x_1..x_3]
Da= x_3^2-x_1^2-x_1*x_2+x_2*x_3
term=x_3^2
pertinentTerm(Da,term)
positivityCheck(o36#0,o36#1)
break
hornTermDa(sub(Da,DualCoord),sub(termsDa#0,DualCoord))
positivityCheck(B,{-1,-1,-1})
exponentiateVector({1,2},{1,3})
positivityCheck(o28#0,o28#1)
positivityCheck(o72#0,flatten entries o72#1)
vec=o28#1
exponentsVec={1,2,3,4,5,6}
apply(#vec,i->(sub(vec#i,RR))^(exponentsVec#i))

break
-- TO DO: Reproduce the results,
-- is the function pertinentDiscriminant really taking care of the positivity conditions?
-- Maybe we need to check again.
-- Make sure we deal correctly with the rings.
-- DualCoords is hardwired so it is important to work in the correct ring.
