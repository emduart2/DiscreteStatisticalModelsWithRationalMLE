
restart
-- CODE INSTRUCTIONS:
-- Run the following command in Mcaulay2 to load all the functions. Change the path according to where you have placed
-- the functionsMLE.m2 file

load "~/Documents/gitHubRepos/DiscreteStatisticalModelsWithRationalMLE/functionsMLE.m2"

-- After you load the files above you can start reproducing the examples considered in https://arxiv.org/abs/1903.06110.
-- All of these examples are numbered according to their number in https://arxiv.org/abs/1903.0611.
-- To reproduce each of the examples simply run each of the lines in Macaulay2.

-- List of functions = {hUniform,toFraction,psiCoord,checkFriendly,checkPair,implicitEquations,pertinentTerm,
-- pertinentDiscriminant,toMonomial,hornTermDa,toHmonomials,pertinentB}
-- Each function in functionsMLE.m2 has a short description of the computations it performs.


------------------------
-- EXAMPLE 5
------------------------
-- A-Discriminant of a cubic
A         =  matrix{{0,1,2,3}}
DualCoord = QQ[x_1,x_2,x_3,x_4]
Da        =  sub(aDiscriminantA A, DualCoord) -- Computes the A-discriminant
resultsDa = pertinentDiscriminant(sub(Da,DualCoord)) -- Checks what terms of the A-discriminat give a satistical model
term=-27*x_1^2*x_4^2 -- This term of the discriminant gives a model
resultsDa#4  -- This is the Horn pair, with the marked monomial from Da that gives a statistical model.


----------------------
-- Example 25. 
----------------------
-- For distinct positive integers alpha,beta,gamma with gcd(alpha,beta,gamma)=1, let
-- A = (0,alpha,beta,gamma).
----------------------
-- STEP 1:
----------------------
-- Option 1: Load all discriminants from the saved file "trinomialDiscriminants"
DualCoord=QQ[x_1..x_4]
trinomialDiscrims = value get "trinomialDiscrims";
-- Jump to STEP 2.
------------------------
-- Option 2: Re calculate the A-discriminants.
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
#Amatrices  -- There are 613 elements in Amatrices.
--- 
--- In the next step we compute the discriminant of all the entries in Amatrices.
--- This is long. The user can skip this step and use the saved list trinomialDiscriminants
--- It is important to specify what is the ring DualCoord in which the polynomials
--- in the list trinomialDiscrims live. This is important to
--- call the function pertinentDiscriminant.
DualCoord=QQ[x_1..x_4]
-- Option 1: Compute all discriminants from matrices in the list Amatrics
trinomialDiscrims = time apply(Amatrices,i->sub(aDiscriminantA(matrix{i}),DualCoord));
"trinomialDiscrims" << toString trinomialDiscrims << endl << close -- saved output to file
#trinomialDiscrims

----------------------
-- STEP 2:
----------------------
-- Need to make sure to load the correct DualCoord ring so that the code works.
DualCoord=QQ[x_1..x_4]
-- This code applys the function pertinentDiscriminant to each entry in trinomialDiscriminants. The output of pertinentDiscriminant
-- is a list that indicates which terms in the tested discriminat give statistical models
testTrinomialDiscrims= for j from 0 to #trinomialDiscrims-1 list (pertinentDiscriminant(sub(trinomialDiscrims#j,DualCoord)));
#(flatten testTrinomialDiscrims) -- This the total number of tested pairs
rawResults = delete("Not pertinent",flatten testTrinomialDiscrims);
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
for i from 0 to 10 list statisticalModels#i

----------------------
-- Example 26. 
----------------------
-- For distinct positive integers alpha,beta,gamma, epsilon  with gcd(alpha,beta)= gcd(gamma,epsilon)=1, let
-- A = {{0,alpha,beta,0,gamma,epsilon},{0,0,0,1,1,1}}.
----------------------
-- STEP 1:
----------------------
-- Option 1: Load all discriminants from the saved file "resultantDiscriminants"
DualCoord=QQ[x_1..x_6]
resultantDiscrims = value get "resultantDiscriminants";
-- Jump to STEP 2.
-- Option 2: Compute the A-discriminants of al matrices A.
-- First we list all pairs and then compute the discriminants
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
resultantDiscrims =  time for k from 0  to #Vpairs-1 list(sub(aDiscriminantA((matrix Vpairs#k)||matrix{{0,0,0,1,1,1}}),DualCoord)); -- 1016.41 seconds to compute
#resultantDiscrims-1
"resultantDiscriminants" << toString resultantDiscrims << endl << close -- saved output to file
----------------------
-- STEP 2:
----------------------

DualCoord = QQ[x_1..x_6] -- Make sure the ring is the correct one.
testResultantDiscrims = time for j from 0 to (#resultantDiscrims-1) list (pertinentDiscriminant(sub(resultantDiscrims#j,DualCoord))); --used 7.5 seconds
rawResults = delete("Not pertinent",flatten testResultantDiscrims);
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
for i from 0 to 10 list statisticalModels#i -- List the first ten statistical models
---------------------------------------------------------------------------------------
----------------------
-- Example 27. 
----------------------

DualCoord=QQ[x_1..x_4]
DA = x_1+x_2+x_3+x_4
---------------------------------------
--- Code for postive binomial multiples
---------------------------------------
a=1
deg=9
binomialPlusMultiples={};
--- STEP 1:
--- List all the binomial multiples (x^a + x^b)*DA
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
--- STEP 2:
--- Test which ones are statistical models
testBinomialPlusMultiples = for j from 0 to #binomialPlusMultiples -1 list pertinentDiscriminant(binomialPlusMultiples#j);
#(flatten testBinomialPlusMultiples)
rawResults = delete("Not pertinent",flatten testBinomialPlusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
statisticalModels
"ModelsBinomialPlusMultiples" << toString statisticalModels << endl << close -- saved output to file
-----------------------------------------
--- Code for negative binomial multiples
---------------------------------------
a=1
deg=9
binomialMinusMultiples={};
--- STEP 1:
--- List all the binomial multiples (x^a - x^b)*DA
while a < deg do(
bzz= flatten entries super basis(a,DualCoord);
binomials = subsets(bzz,2);
currentBinomialMultiple = for j from 0 to #binomials-1 list( if(gcd(binomials#j)!=1) then continue ; DA*(sum {(binomials#j)#0,-1*(binomials#j)#1}));
print(#currentBinomialMultiple);
binomialMinusMultiples=binomialMinusMultiples|currentBinomialMultiple;
a=a+1;
)
#binomialMinusMultiples
--- STEP 2:
--- Test which ones are statistical models
testBinomialMinusMultiples = for j from 0 to #binomialPlusMultiples -1 list pertinentDiscriminant(binomialMinusMultiples#j);
#(flatten testBinomialMinusMultiples)
rawResults = delete("Not pertinent",flatten testBinomialMinusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
statisticalModels
"ModelsBinomialMinusMultiples" << toString statisticalModels << endl << close -- saved output to file
---------------------------------------------
--- Code for positive trinomial multiples
---------------------------------------------
a=1
deg=4
trinomialPlusMultiples={};
--- STEP 1:
--- List all the trinomial multiples (x^a + x^b + x^c)*DA
while a < deg do(
bzz= flatten entries super basis(a,DualCoord);
trinomials = subsets(bzz,3);
currentBinomialMultiple = for j from 0 to #trinomials-1 list( if(gcd(trinomials#j)!=1) then continue ; DA * (sum trinomials#j));
print(#currentBinomialMultiple);
trinomialPlusMultiples=trinomialPlusMultiples|currentBinomialMultiple;
a=a+1;
)
#trinomialPlusMultiples
--- STEP 2:
--- Test which ones are statistical models
testTrinomialPlusMultiples = for j from 0 to #trinomialPlusMultiples -1 list pertinentDiscriminant(trinomialPlusMultiples#j);
#(flatten testTrinomialPlusMultiples)
rawResults = delete("Not pertinent",flatten testTrinomialPlusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
statisticalModels
"ModelsTrinomialPlusMultiples" << toString statisticalModels << endl << close -- saved output to file
---------------------------------------------
--- Code for negative trinomial multiples
---------------------------------------------
a=1
deg=4
trinomialMinusMultiples={};
--- STEP 1:
--- List all the trinomial multiples (x^a + x^b + x^c)*DA
while a < deg do(
bzz= flatten entries super basis(a,DualCoord);
trinomials = subsets(bzz,3);
currentBinomialMultiple = for j from 0 to #trinomials-1 list( if(gcd(trinomials#j)!=1) then continue ; DA * (sum {(trinomials#j)#0,(trinomials#j)#1,-1*(trinomials#j)#2}));
print(#currentBinomialMultiple);
trinomialMinusMultiples=trinomialMinusMultiples|currentBinomialMultiple;
a=a+1;
)
#trinomialMinusMultiples
--- STEP 2:
--- Test which ones are statistical models
testTrinomialMinusMultiples = for j from 0 to #trinomialMinusMultiples -1 list pertinentDiscriminant(trinomialMinusMultiples#j);
#(flatten testTrinomialMinusMultiples)
rawResults = delete("Not pertinent",flatten testTrinomialMinusMultiples);
#rawResults
statisticalModels= for j from 0 to #rawResults-1 list (if(last rawResults#j== false) then continue;rawResults#j);
#statisticalModels
"ModelsTrinomialMinusMultiples" << toString statisticalModels << endl << close -- saved output to file
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
