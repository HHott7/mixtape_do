************************* Mixtape **************************

******** Chapter 1 
** Simulando um banco de dados para aplicar o estimador OLS (ols.do)
set seed 1 
clear 
set obs 10000 
gen x = rnormal() 
gen u  = rnormal() 
gen y  = 5.5*x + 12*u 
reg y x 
predict yhat1 
gen yhat2 = -0.0750109  + 5.598296*x // Compare yhat1 and yhat2
sum yhat* 
predict uhat1, residual 
gen uhat2=y-yhat2 
sum uhat* 
twoway (lfit y x, lcolor(black) lwidth(medium)) (scatter y x, mcolor(black) ///
msize(tiny) msymbol(point)), title(OLS Regression Line) 
rvfplot, yline(0) 

** ols_2.do
clear 
set seed 1234
set obs 10
gen x = 9*rnormal() 
gen u  = 36*rnormal() 
gen y  = 3 + 2*x + u
reg y x
predict yhat
predict residuals, residual
su residuals
list
collapse (sum) x u y yhat residuals //note que a soma dos residuos e zero. Logo y=y_hat na media (na media o erro e zero)
list
// Note que covariancia (somatorio x*residuals) e zero. Alem disso, y_hat e residuals sao ortogonais (
gen ux = x*residuals
gen uy = yhat*residuals

** ols_3.do
clear all 
program define ols, rclass 
version 14.2 
syntax [, obs(integer 1) mu(real 0) sigma(real 1) ] 

    clear 
    drop _all 
    set obs 10000 
    gen x = 9*rnormal()  
    gen u  = 36*rnormal()  
    gen y  = 3 + 2*x + u 
    reg y x 
    end 
// Vamos pega uma amostra com 10 mil observacoes, com x e u seguindo uma normal (cada vez que rodarmos os valores 
//serao diferentes pois nao especificamos uma seed). A ideia e que na media as estimativas serao igual ao beta populacional (2)
simulate beta=_b[x], reps(1000): ols
su 
hist beta

** reganat.do -reganat- is a user-created package by @Filoso2013. Utilizado para demonstrar a importancia do teorema FWL
ssc install reganat, replace
sysuse auto.dta, replace
regress price length //bivariada
regress price length weight headroom mpg //multivariada
reganat price length weight headroom mpg, dis(length) biline

///// Cluster

// Cluster I (estimativas e intervalo de confianca com erros nao clusterizados)
** cluster1.do
clear all
set seed 20140
* Set the number of simulations
local n_sims  = 1000
set obs `n_sims'

* Create the variables that will contain the results of each simulation
generate beta_0 = .     //estimativa de beta0
generate beta_0_l = .   //intervalo de confianca de beta0 (por baixo)
generate beta_0_u = .   //intervalo de confianca de beta0 (por cima)
generate beta_1 = .
generate beta_1_l = .
generate beta_1_u = .


* Provide the true population parameters
local beta_0_true = 0.4
local beta_1_true = 0
local rho = 0.5

* Run the linear regression 1000 times and save the parameters beta_0 and beta_1
quietly {
    forvalues i = 1(1) `n_sims' {
        preserve  
        clear
        set obs 100
        generate x = rnormal(0,1)
        generate e = rnormal(0, sqrt(1 - `rho'))
        generate y = `beta_0_true' + `beta_1_true'*x + e
        regress y x
        local b0 = _b[_cons]
        local b1 = _b[x]
        local df = e(df_r)
        local critical_value = invt(`df', 0.975)
        restore
        replace beta_0 = `b0' in `i'
        replace beta_0_l = beta_0 - `critical_value'*_se[_cons] 
        replace beta_0_u = beta_0 + `critical_value'*_se[_cons] 
        replace beta_1 = `b1' in `i'
        replace beta_1_l = beta_1 - `critical_value'*_se[x] 
        replace beta_1_u = beta_1 + `critical_value'*_se[x] 
        
    }
}
gen false = (beta_1_l > 0 ) //Situacao em que o intervalo de confianca nao contem o parametro populacional (beta1=0)
replace false = 2 if beta_1_u < 0 //Situacao em que o intervalo de confianca nao contem o parametro populacional (beta1=0)
replace false = 3 if false == 0 //Situacao em que o intervalo de confianca contem o parametro populacional (beta1=0)
tab false //95% das vezes o intervalo de confianca contem o beta1 (incorretamente rejeita que beta1=0 apenas 5% das vezes)

* Plot the parameter estimate
hist beta_1, frequency addplot(pci 0 0 100 0) title("Least squares estimates of non-clustered data") subtitle(" Monte Carlo simulation of the slope") legend(label(1 "Distribution of least squares estimates") label(2 "True population parameter")) xtitle("Parameter estimate") 

sort beta_1
gen int sim_ID = _n
gen beta_1_True = 0
* Plot of the Confidence Interval
twoway rcap beta_1_l beta_1_u sim_ID if beta_1_l > 0 | beta_1_u < 0  , horizontal lcolor(pink) || || ///
rcap beta_1_l beta_1_u sim_ID if beta_1_l < 0 & beta_1_u > 0 , horizontal ysc(r(0)) || || ///
connected sim_ID beta_1 || || ///
line sim_ID beta_1_True, lpattern(dash) lcolor(black) lwidth(1) ///  
title("Least squares estimates of non-clustered data") subtitle(" 95% Confidence interval of the slope") ///
legend(label(1 "Missed") label(2 "Hit") label(3 "OLS estimates") label(4 "True population parameter")) xtitle("Parameter estimates") ///
ytitle("Simulation")

// Cluster II (estimativas e intervalo de confianca com erros clusterizados - ignorando o problema)
* cluster2.do
clear all
set seed 20140
local n_sims = 1000
set obs `n_sims'

* Create the variables that will contain the results of each simulation
generate beta_0 = .
generate beta_0_l = .
generate beta_0_u = .
generate beta_1 = .
generate beta_1_l = .
generate beta_1_u = .


* Provide the true population parameters
local beta_0_true = 0.4
local beta_1_true = 0
local rho = 0.5

* Simulate a linear regression. Clustered data (x and e are clustered)

quietly {
forvalues i = 1(1) `n_sims' {
    preserve
    clear
    set obs 50
    
    * Generate cluster level data: clustered x and e
    generate int cluster_ID = _n
    generate x_cluster = rnormal(0,1)
    generate e_cluster = rnormal(0, sqrt(`rho'))
    expand 20
    bysort cluster_ID : gen int ind_in_clusterID = _n

    * Generate individual level data
    generate x_individual = rnormal(0,1)
    generate e_individual = rnormal(0,sqrt(1 - `rho'))

    * Generate x and e
    generate x = x_individual + x_cluster
    generate e = e_individual + e_cluster
    generate y = `beta_0_true' + `beta_1_true'*x + e
    
* Least Squares Estimates
    regress y x
    local b0 = _b[_cons]
    local b1 = _b[x]
    local df = e(df_r)
    local critical_value = invt(`df', 0.975)
    * Save the results
    restore
    replace beta_0 = `b0' in `i'
    replace beta_0_l = beta_0 - `critical_value'*_se[_cons]
    replace beta_0_u = beta_0 + `critical_value'*_se[_cons]
    replace beta_1 = `b1' in `i'
    replace beta_1_l = beta_1 - `critical_value'*_se[x]
    replace beta_1_u = beta_1 + `critical_value'*_se[x]
}
}

gen false = (beta_1_l > 0 ) //Situacao em que o intervalo de confianca nao contem o parametro populacional (beta1=0)
replace false = 2 if beta_1_u < 0 //Situacao em que o intervalo de confianca nao contem o parametro populacional (beta1=0)
replace false = 3 if false == 0 //Situacao em que o intervalo de confianca contem o parametro populacional (beta1=0)
tab false //Com o erro custerizado o intervalo de confianca ja nao contem o verdadeiro beta1 em 95% das vezes (apenas 58%) - erro tipo1 aumentou muito

* Plot the parameter estimate
hist beta_1, frequency addplot(pci 0 0 100 0) title("Least squares estimates of clustered Data") subtitle(" Monte Carlo simulation of the slope") legend(label(1 "Distribution of least squares estimates") label(2 "True population parameter")) xtitle("Parameter estimate")
//A estimativa do parametro nao e afetada, ja que estamos mudando apenas a inferencia

* cluster3.do
sort beta_1
gen int sim_ID = _n
gen beta_1_True = 0

* Plot of the Confidence Interval
twoway rcap beta_1_l beta_1_u sim_ID if beta_1_l > 0 | beta_1_u < 0  , horizontal lcolor(pink) || || ///
rcap beta_1_l beta_1_u sim_ID if beta_1_l < 0 & beta_1_u > 0 , horizontal ysc(r(0)) || || ///
connected sim_ID beta_1 || || ///
line sim_ID beta_1_True, lpattern(dash) lcolor(black) lwidth(1) ///  
title("Least squares estimates of clustered data") subtitle(" 95% Confidence interval of the slope") ///
legend(label(1 "Missed") label(2 "Hit") label(3 "OLS estimates") label(4 "True population parameter")) xtitle("Parameter estimates") ///
ytitle("Simulation")

** Veja que as vezes em que nao encontramos o verdadeiro beta1 aumentam significativamente (era 5% agora 42% das vezes, como vimos no tab false)
** Ou seja, incorretamente rejeitamos a hipotese de que beta1=0 muitas vezes quando os dados sao clusterizados e ignoramos isso

// Cluster III (estimativas e intervalo de confianca com erros clusterizados - lidando com o problema)
* cluster4.do
* Robust Estimates
clear all
local n_sims = 1000
set obs `n_sims'

* Create the variables that will contain the results of each simulation
generate beta_0_robust = .
generate beta_0_l_robust = .
generate beta_0_u_robust = .
generate beta_1_robust = .
generate beta_1_l_robust = .
generate beta_1_u_robust = .

* Provide the true population parameters
local beta_0_true = 0.4
local beta_1_true = 0
local rho = 0.5

quietly {
forvalues i = 1(1) `n_sims' {
    preserve
    clear
    set obs 50
    
    * Generate cluster level data: clustered x and e
    generate int cluster_ID = _n
    generate x_cluster = rnormal(0,1)
    generate e_cluster = rnormal(0, sqrt(`rho'))
    expand 20
    bysort cluster_ID : gen int ind_in_clusterID = _n

    * Generate individual level data
    generate x_individual = rnormal(0,1)
    generate e_individual = rnormal(0,sqrt(1 - `rho'))

    * Generate x and e
    generate x = x_individual + x_cluster
    generate e = e_individual + e_cluster
    generate y = `beta_0_true' + `beta_1_true'*x + e
    regress y x, cl(cluster_ID) //agora regride considerando a clusterizacao
    local b0_robust = _b[_cons]
    local b1_robust = _b[x]
    local df = e(df_r)
    local critical_value = invt(`df', 0.975)
    * Save the results
    restore
    replace beta_0_robust = `b0_robust' in `i'
    replace beta_0_l_robust = beta_0_robust - `critical_value'*_se[_cons]
    replace beta_0_u_robust = beta_0_robust + `critical_value'*_se[_cons]
    replace beta_1_robust = `b1_robust' in `i'
    replace beta_1_l_robust = beta_1_robust - `critical_value'*_se[x]
    replace beta_1_u_robust = beta_1_robust + `critical_value'*_se[x]

}
}

* Plot the histogram of the parameters estimates of the robust least squares
gen false = (beta_1_l_robust > 0 )
replace false = 2 if beta_1_u_robust < 0
replace false = 3 if false == 0
tab false //Mesmo com clusterizacao agora o intervalo de confianca contem o verdadeiro parametro 98% das vezes (ou seja, o erro tipo 1 caiu muito)

* Plot the parameter estimate
hist beta_1_robust, frequency addplot(pci 0 0 110 0) title("Robust least squares estimates of clustered data") subtitle(" Monte Carlo simulation of the slope") legend(label(1 "Distribution of robust least squares estimates") label(2 "True population parameter")) xtitle("Parameter estimate")

sort beta_1_robust
gen int sim_ID = _n
gen beta_1_True = 0

* Plot of the Confidence Interval
twoway rcap beta_1_l_robust beta_1_u_robust sim_ID if beta_1_l_robust > 0 | beta_1_u_robust < 0, horizontal lcolor(pink) || || rcap beta_1_l_robust beta_1_u_robust sim_ID if beta_1_l_robust < 0 & beta_1_u_robust > 0 , horizontal ysc(r(0)) || || connected sim_ID beta_1_robust || || line sim_ID beta_1_True, lpattern(dash) lcolor(black) lwidth(1) title("Robust least squares estimates of clustered data") subtitle(" 95% Confidence interval of the slope") legend(label(1 "Missed") label(2 "Hit") label(3 "Robust estimates") label(4 "True population parameter")) xtitle("Parameter estimates") ytitle("Simulation")
