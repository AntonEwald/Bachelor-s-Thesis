Kandidat Skrivandet
================
Anton Holm
2020-02-24

## Abstrakt

## Introduction

## Background

At the Swedish Museum of Natural History, the Department of
Environmental Research and Monitoring in a joint effort with other
departments conducts statistical research of environmental toxicants as
part of the National Swedish Contaminant Programme in marine biota. One
of the programs conducted regards analysing long term time trends of
several toxins in Swedish waters and to estimate the rate of change. The
models used to analyse such time trends are at the moment surprisingly
elemental and disregards much of the data collected. One of the more
common, but nonetheless crucial oversights, concerns building models and
drawing conclusions from fabricated data due to data being censored.

## Data

The report from Bignert at al (2017) explains much of the data sampling.
The data comes from several sampling areas regarded as locally
uncontaminated. Several species of fish, as well as guillemot eggs and
blue mussels, are collected from different sampling areas each year.
When collected, a constant number of 10-12 speciemens independent of
each other are analysed for a large number of toxins. For some species,
the analysis is done for pooled samples containing a number of
speciemens in each pool. To reduce the between-years variation, each
sampling area tries to analyse speciemens of the same sex and age.
However, the variation can not be reduced to zero and other parameters
effects the variation such as fat content and local discharges as an
example. The concentration between each fish will also contain noise,
hence the data sampled will have variation between years as well as
within years.

As a result of test equipments not being able to detect small enough
quantities of toxins, a portion of the data is reported as *below the
limit of quantification (LOQ)*. This portion of the data is reported as
the LOQ divided by the square root of 2.

Due to biological properies such as size and fat tissues being able to
effect the concentration of toxins and these attributes being effected
by sampling site, this thesis will analyse sampling areas individually.

Bignert at al (2017) uses log-linear regression analysis, hence the data
is assumed to follow a log-linear distribution.

## Common Errors

One of the most common error being made when analysing censored data is
fabricating. The analysists simply substitute the non-detects with a
fraction (often one half) of the quantitive- or detetection limit. A
simulation were made by Helsel (2006) showcasing that this method
produces lousy estimates of statistics and have the potential to not
only overlook patterns in the data, but also impose it’s own fabricated
patterns. This could result in a goverment investing millions to clean a
lake of toxins after a report displaying an increase in concentrations
of a certain metal in fish when in fact, there were no such pattern to
begin with. The reverse is even more terrifying, obtaining a report
showing no significant increase in concentration, when indeed the
concentration of said metal have been increasing for years. Causes of an
increase in concentration have been missed, remediations goes undone and
the health of humans and the ecosystem is unnecessarily endangered.
There are plenty more mistakes commonly being made when handling
censored data including misinterpreting an improvement in measuring
technique for a decrease in non-detects. However, this will not be
discussed in detail in this thesis.

# Theory

When working with censored data, the non-detects can’t be looked at as
having a specific value. Instead, a combination of the information of
the proportion of non-detects with the numerical values of the
uncensored observations gives a better understanding of the data.
Assuming a distribution for the data above and below the reported limit
in combination with the above mentioned information gives a foundation
to work with maximum likelihood estimates (MLE). In a study of Chung
(1990) regarding regression analysis of geochemical data with
non-detects, it was shown that MLE gave a much better estimation for the
true value of the slope coefficient than any of the substitution values
(0, 0.1, …, 1 times the detection limit). Regression analysis for
censored data is being used in many fields, including but not limited
to, medical statistics as used by Lee and Go (1997) and in economics
where Chay and Honore (1998) used MLE regression on right-censored data
to model incomes. However, for left-censored data where the residuals is
assumed to follow a normal distribution, the MLE regression is sometimes
mentioned as Tobit analysis after the famous economist James Tobin. For
the particular data from the Museum of Natural History, the use of Tobit
regression models can serve useful to handle the censoring while the use
of a Linear Mixed-Effect Model (LMM) will deal with the fact that data
contains variation both within and between years.

## CDF of a linear regression model

Consider a normal simple linear regression model \[
y_i = x_i \beta + \epsilon_i, \; \epsilon_i \sim N(0,\sigma^2)
\]

were \(y_i\) is the response variable, \(x_i\) the explanatory variable,
\(\beta\) an effect parameter and \(\epsilon_i\) the error term. It’s
then easy to find the cumulative distribution function (CDF) for this
model.

\[
F(y_i) = P(x_i\beta+\epsilon_i\leq y_i) = P(\frac{\epsilon_i}{\sigma}\leq \frac{1}{\sigma}(y_i-x_i\beta)) = \Phi[\frac{1}{\sigma}(y_i-x_i\beta)]
\]

where \(\Phi(\cdot)\) is the CDF for a standard normal variable. The
probability density function (PDF) is further given by
\(f(y_i)=\frac{dF(y_i)}{dy_i}\).

## Linear mixed-effects model

\(\mathbf{**}\)Can aggregate data. Take mean of each group =\> the avg
data points are now independent: Less noise but disregard a lot of data
Can do regression on each group =\> a lot of noise but takes all data
LMM somewhere in between\(\mathbf{**}\)

Mixed models are an extension of normal models where random effects are
integrated. A linear mixed model is an extension of mixed models where
both the fixed and random effects take place linearly in the model. The
random effects can be observed as additional error terms in the model.
Following the notation of Pinheiro and Bates (2000) the linear mixed
model for a single level of grouping, as described by Laird and Ware
(1982), can be expressed as

\[
\mathbf{y_i} = \mathbf{X_i}\mathbf{\beta} + \mathbf{Z_i}\mathbf{b_i} + \mathbf{\epsilon_i}
\]

for \(i = 1,...,M\). Here, \(\mathbf{y_i}\) is the \(n_i\) dimension
respons vector for group i, \(\mathbf{\beta}\) the \(p\) dimensional
vector of fixed-effect parameters, \(\mathbf{b_i}\) the \(q\)
dimensional vector of random-effects, \(\mathbf{X_i}\) a matrix with
covariates of size \(n_i\) x \(p\), \(\mathbf{Z_i}\) a design matrix of
size \(n_i\) x \(q\) linking \(\mathbf{b_i}\) to \(\mathbf{y_i}\) and
\(\mathbf{\epsilon_i}\) an \(n_i\) dimension vector of error terms
within group i with \(\mathbf{b_i}\sim N(0,\Sigma)\), \(\Sigma\) being
the symmetrical, positive semi-definite \(n_i\) x \(n_i\) dimension
covariance matrix and \(\mathbf{\epsilon_i}\sim N(0,\sigma^2I)\), \(I\)
being the \(n_i\) dimension vector of ones.

## Tobit Model

The Tobit model is characterized by the latent regression equation \[
y_i^* = \mathbf{x_i}\cdot\mathbf{\beta} + \epsilon_i, \; \epsilon_i \sim N(0, \sigma^2)
\]

where \(y_i^*\) is the laten dependent variable, \(\mathbf{x_i}\) is a
vector of covariates, \(\mathbf{\beta}\) a vector of effect parameters
and \(\epsilon_i\) is the error term. Given this, the observed dependent
variable can be specified as:

\[
\begin{cases}
y_i = y_i^*, & y_i^* > y_L \\
y_i = y_L, & otherwise
\end{cases}
\]

with \(y_L\) being the reporting limit. This leads us to the PDF of the
Tobit model:

\[
f(y_i|\mathbf{x_i}) = \begin{cases}
f(y_i|\mathbf{x_i}) = 0, & y_i<y_L\\
f(y_L|\mathbf{x_i}) = P(y_i^* \leq y_L|\mathbf{x_i}), & y_i=y_L\\
 f(y_i|\mathbf{x_i})=f(y_i^*|\mathbf{x_i}), & y_i>y_L
\end{cases}
\]

Using the same method as for a normal simple linear regression model, we
further deduce

\[
f(y_i|x_i)= \begin{cases}
0, & y_i<y_L \\
\Phi\bigg{(}\frac{y_L-\mathbf{x_i}\cdot\mathbf{\beta}}{\sigma}\bigg{)}, & y_i=y_L \\
\frac{1}{\sigma}\phi\bigg{(}\frac{y_i-\mathbf{x_i} \cdot \mathbf{\beta}}{\sigma}\bigg{)}, & y_i>y_L
\end{cases}
\]

where \(\phi(\cdot)\) is the PDF of a standard normal distribution.
Hence, the likelihood function for the Tobit model is:

\[
L = \underset{y_i=y_L}{\Pi} \Phi\bigg{(}\frac{y_L-\mathbf{x_i}\cdot\mathbf{\beta}}{\sigma}\bigg{)} \cdot \underset{y_i>y_L}{\Pi}\frac{1}{\sigma}\phi\bigg{(}\frac{y_i-\mathbf{x_i}\cdot\mathbf{\beta}}{\sigma}\bigg{)}
\]

## Case of the museum (Multivariate Normal-distribution)

Now, in the case of the analysis conducted by the Swedish Museum of
Natural History, a Linear Mixed Tobit Model could be implemented.
Regarding each year as a seperate group \(t\) having \(n_t\) speciemens.
The between-year variance is the same for each speciemen in the same
group while the within-year variance is the same for every speciemen
through each year.

Hence, the model is

\[
\log(\mathbf{y_t}) = \mathbf{x_t} \cdot \mathbf{\beta} + \mathbf{z} \cdot \mathbf{e_t} + \mathbf{\epsilon}
\]

where \(\mathbf{y_t}\) is the \(n_t\) dimension response vector
containing the measured concentration of a certain toxin,
\(\mathbf{x_t}\) a matrix of dimension \(n_t\) x \(2\) having a column
of ones for the intercept and a column of the year of sampling,
\(\mathbf{\beta}\) the 2 dimensional vector of fixed effect parameters
including the intercept, \(\mathbf{z}\) an \(n_t\) dimensional row
vector of ones, \(\mathbf{e_t}\) an \(n_t\) dimensional vector of the
random effect \(e_t\) and \(\mathbf{\epsilon}\) the \(n_t\) dimensional
vector with the within-years variance for each speciemen
\(\epsilon_i, i=1,2,...,n_i\). Further more, since
\(e_t \sim N(0,\sigma_t^2)\) and \(\epsilon \sim N(0,\delta^2)\), the
distribution of \(\log(\mathbf{y_t})\) follows

\[
\log(\mathbf{y_t}) \sim N_{n_t}(\mathbf{x_t}\cdot \mathbf{\beta}, \mathbf{\Sigma})
\]

with \(\mathbf{\Sigma} = (a_{ij})\in \mathbb{R}^{n_t \text{x} n_t}\) the
covariance matrix where
\((a_{ij}) = Cov (e + \epsilon_i, e +\epsilon_j)\). Further calculations
of the covariance gives

\[
Cov(e + \epsilon_i, e + \epsilon_j) = E[(e+\epsilon_i)(e+\epsilon_j)] - E[e+\epsilon_i]E[e+\epsilon_j]= E[\epsilon^2] = \delta^2
\]

for all \(i,j\) such that \(i\neq j\) since \(E[e]=E[\epsilon_k]=0\) for
all \(k\). In addition,
\((a_{ij}) = Var(e +\epsilon_i) = \sigma^2 + \delta^2\) when \(i=j\).

Following the method above used to derive the CDF of a linear regression
model, the CDF of the model in question can also be derived. First of
all, the fact that observations can be censored must be taken into
consideration. This is done by partioning the data into censored and
non-censored components

\[
\mathbf{y_t}=
    \begin{bmatrix}
           \mathbf{y_t^o} \\
           \mathbf{y_t^c} \end{bmatrix}
           \mathbf{x_t} =  \begin{bmatrix}
           \mathbf{x_t^o} \\
           \mathbf{x_t^c} \end{bmatrix}
         \mathbf{\Sigma_t} = \begin{bmatrix}
           \mathbf{\Sigma_t^{oo}} & \mathbf{\Sigma_t^{oc}}\\
           \mathbf{\Sigma_t^{oc^{T}}} & \mathbf{\Sigma_t^{cc}}
         \end{bmatrix}
\]

where \(\mathbf{y_t^o}\) is the \(n_t^o\) vector of all the observed,
non-censored values and \(\mathbf{y_t^c}\) the \(n_t^c\) vector of all
censored observations, the same following for \(\mathbf{x_t}\) being
partioned into a \(n_t^o \, \text{x} \, 2\) matrix and a
\(n_t^c \, \text{x} \, 2\) matrix while \(\mathbf{\Sigma_t^{oo}}\) and
\(\mathbf{\Sigma_t^{cc}}\) is the matrix of variances and covariances
between all observed values and censored values respectively and
\(\mathbf{\Sigma_t^{oc}}=\mathbf{\Sigma_t^{co^{T}}}\) being the matrix
of covariances between non-censored and censored observations. It
follows that \(\mathbf{y_t^o}\) has a multivariate normal distribution
with PDF \(f_{\mathbf{y_i}^o}\). Using the properties of the
multivariate normal distribution, following Eaton (1983), the
conditional distribution of \(y_t^c|y_t^o\) is also multivariate
normally distributed with mean and variance as follows

\[
\mu_t^{c|o} = \mathbf{x}_t^c\mathbf{\beta} + \mathbf{\Sigma_t^{co}}\mathbf{\Sigma_t^{{oo}^{-1}}}(\mathbf{y_t^o}-\mathbf{x_t^o\beta}), \;\;\;\; \mathbf{\Sigma_t^{c|o}} = \mathbf{\Sigma_t^{cc}}-\mathbf{\Sigma_t^{co}}\mathbf{\Sigma_t^{{oo}^{-1}}}\mathbf{\Sigma_t^{co^{T}}}
\]

here \(\Sigma_t^{oo^{-1}}\) is the inverse of \(\Sigma_t^{oo}\). Denote
\(\phi_t^{c|o}(\cdot)\) as the PDF of the conditional distribution
function of \(y_t^c\) given \(y_t^o\) and \(\mathbf{c_t}\) the \(n_t^c\)
vector where \(c_{tj}\) is the censoring threshold for the \(j^{th}\)
censored outcome. Now, since all \(\mathbf{y_t}\) are independent, using
the methods of previous sections, the likelihood function can be written
as

\[
L(\mathbf{\beta}) = \underset{t}{\Pi} f_{\mathbf{y_t}^o}(\mathbf{y_t}^o|\mathbf{\beta})\cdot \phi_t^{c|o}(\mathbf{c_t}|\mathbf{\beta})
\]

which given the PDF of a multivariate normal distributed variable gives
\[
\begin{aligned}
& \underset{t}{\Pi} \frac{1}{\sqrt{(2\pi)^{n_t^o}|\mathbf{\Sigma_t}^{oo}|}}\cdot exp\bigg{\{}-\frac{1}{2}(\mathbf{y_t}^o-\mathbf{x_t}^o\mathbf{\beta})^T\mathbf{\Sigma}_t^{oo^{-1}}(\mathbf{y_t}^o-\mathbf{x_t}^o\mathbf{\beta})\bigg{\}} \cdot \\
& \int_{-\infty}^{n_{t1}}\int_{-\infty}^{n_{t2}}\cdots\int_{-\infty}^{n_{tn_t^c}} \frac{1}{\sqrt{(2\pi)^{n_t^c}|\mathbf{\Sigma_t}^{c|o}|}}\cdot exp \bigg{\{}-\frac{1}{2}(\mathbf{z}-\mathbf{\mu^{c|o}})^T\mathbf{\Sigma}_t^{c|o^{-1}}(\mathbf{z}-\mathbf{\mu^{c|o}})\bigg{\}}
\end{aligned}
\]

## LMEC

#### EM-Algorithm

## Simulation

We want a slope representing a 1% yearly increase. Our model is
\(Y=e^{\beta_1 X + \epsilon}\) so when \(X\) goes to \(X+1\) we want
\(Y\) to go to \(Y\cdot 1.01\) hence
\(Y(x+1)=e^{\beta_1 (X + 1) + \epsilon} = e^{\beta_1 X + \epsilon}e^{\beta_1}=Y\cdot e^{\beta_1}\)
hence \(e^\beta_1=1.01\) so \(\beta_1=log(1.01)\) so in the log scale,
our slope is log(log(1.01)) : Probobly wrong

## Result

## Conclusion

## References

1)  Helsel, D.R., 2006, Fabricating data: how substituting values for
    censored observations can ruin results, and what can be done about
    it. Chemosphere 65, pp. 2434–2439, doi:
    <https://doi.org/10.1016/j.chemosphere.2006.04.051>

2)  Chung, C.F., 1990, Regression analysis of geochemical data with
    observations below detection limits, in G. Gaal and D.F.
    Merriam,eds., Computer Applications in Resource Estimation.
    Pergammon Press, New York, pp. 421–433, doi:
    <https://doi.org/10.1016/B978-0-08-037245-7.50032-9>

3)  Lee, T.L and Go, O.T, 1997, Survival Analysis in Public Health
    Research, vol.18, pp. 105-134, doi:
    <https://doi.org/10.1146/annurev.publhealth.18.1.105>

4)  Chay, K.Y. and Honore, B.E. , 1998, Estimation of censored
    semiparametric regression models: an application to changes in
    Black–White earnings inequality during the 1960s. Journal of Human
    Resources Vol.33, pp. 4–38, doi: 10.2307/146313

5)  Pinheiro, J.C and Bates, D.M, (2000), Mixed-Effects Models in S and
    S-PLUS (1. ed.), New York: Springer

6)  Laird, N. M. and Ware, J. H. (1982). Random-effects models for
    longitudinal data, Biometrics 38: 963–974.

7)  Bignert, A., Danielsson, S., Faxneld, S., Ek, C., Nyberg, E. (2017).
    Comments Concerning the National Swedish Contaminant Monitoring
    Programme in Marine Biota, 2017, 4:2017, Swedish Museum of Natural
    History, Stockholm, Sweden, Retrieved from the website of the Museum
    of Natural Historys:
    <http://nrm.diva-portal.org/smash/get/diva2:1090746/FULLTEXT01.pdf>

8)  Eaton, M. L. (1983). Multivariate Statistics: a Vector Space
    Approach. John Wiley and Sons. pp. 116–117. ISBN 978-0-471-02776-8
