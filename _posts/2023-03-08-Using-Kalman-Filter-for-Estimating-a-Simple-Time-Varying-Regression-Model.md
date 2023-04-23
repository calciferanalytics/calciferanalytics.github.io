---
layout: post
title: "Derivation of the Kalman Filtering Equations for a Time-varying Intercept Simple Linear Regression Model"
date: 2023-04-22
use_math: true
---

Some of the most useful and interesting applications of Kalman filter in economics and finance are in the context of multivariate models. However, the method is best understood when introduced in the context of a simple, univariate model. Therefore, the goal of this post is to introduce the reader to the concept of the Kalman filter using a very simple, univariate regression model. 

Suppose, we have the following simple regression model for the GDP growth:
$$\begin{align}
y_t = \beta_t + e_t \tag{1}
\end{align}$$
where $y_{t}$ is the growth rate of real US GDP and $$e_{t}\sim i.i.d.\text{ } N(0,R). \tag{2}$$ and $\beta_t$ is the time-varying mean GDP growth that follows an AR(1) process:
$$\begin{equation}
\beta _{t}=\mu +F\beta _{t-1}+v_{t} \tag{3}
\end{equation}
$$
where 
$$
v_{t}\sim i.i.d.\text{ } N(0,Q) \tag{4}
$$

The idea of Kalman FIlter is pretty straight-forward and consists of two steps:
1. **Prediction:** At the beginning of time $t$ we may want to form an
optimal predictor of $y_{t}$ based on all the available information up to
time $t-1:$ $y_{t|t-1}$. To do this we need to calculate \ $\beta _{t|t-1}$
2.  **Updating**: Once $y_{t}$ is realized at the end of time $t$, the
prediction error can be calculated: $\eta _{t|t-1}=y_{t}-y_{t|t-1}$. This
prediction error contains new information about $\beta _{t}$ beyond that
contained in $\beta _{t|t-1}$. Thus, after observing $y_{t}$ a more accurate
inference can be made of $\beta _{t}:\beta _{t|t}$. An inference of $\beta_t$ based on information up to time $t$ may be of the following form: $\beta _{t|t}=\beta _{t|t-1}+K_{t}\eta _{t|t-1}$ where $K_t$, the Kalman gain, is the weight assigned to new information about $\beta _{t}$ contained in the
prediction error. To be more specific the basic filter is described in
following two steps:

1) **Prediction** 
$$
\begin{align}
\beta _{t|t-1}=\mu +F\beta _{t-1|t-1} \tag{5}\\
P_{t|t-1}=F^{2}P_{t-1|t-1}+Q  \tag{6}\\
\eta _{t|t-1}=y_{t}-y_{t|t-1}=y_{t}-\beta _{t|t-1} \tag{7}\\
f_{t|t-1}=P_{t|t-1}+R  \tag{8}\\
\end{align}
$$

where 
* $\beta _{t|t-1}=E(\beta _{t}|\psi _{t-1})$ - estimate of $\beta _{t}$ conditional on information up to $t-1$
* $\beta _{t|t}=E(\beta _{t}|\psi _{t})$ - estimate of $\beta _{t}$ conditional on information up to $t$
* $P_{t|t-1}=E(\beta _{t}-\beta _{t|t-1})^{2}$ \ - variance of $\beta _{t}$ conditional on information up to $t-1$
* $P_{t|t}=E(\beta _{t}-\beta _{t|t})^{2}$ \ - variance of $\beta _{t\text{ }}$ conditional on information up to $t$
* $\eta _{t|t-1}=y_{t}-y_{t|t-1}$ - prediction error
* $f_{t|t-1}=E[\eta _{t|t-1}^{2}]$ - conditional variance of prediction error

2. **Updating**
$$ \begin{align}
K_{t}=\frac{P_{t|t-1}}{f_{t|t-1}} \tag{9} \\
\beta _{t|t}=\beta _{t|t-1}+K_{t}\eta _{t|t-1} \tag{10}  \\
P_{t|t}=P_{t|t-1}-KP_{t|t-1} \tag{11} 
\end{align}
$$


Provided that eigenvalues of $F$ are all inside the unit circle, then the
process for $\beta _{t}$ in equation (3) is is covariance-stationary and thus the derivation of equation (5) is straight-forward. Equation (6) can be derived as follows:

$$ \begin{align}
P_{t|t-1} &=E_{t-1}(\beta _{t}-\beta _{t|t-1})^{2}
\\
&=E_{t-1}(\mu +F\beta _{t-1}+v_{t}-\mu -F\beta _{t-1|t-1})^{2}  \\
&=E_{t-1}(F\beta _{t-1}+v_{t}-F\beta _{t-1|t-1})^{2}   \\
&=E_{t-1}(F(\beta _{t-1}-\beta _{t-1|t-1})+v_{t})^{2}  \\
&=F^{2}E_{t-1}(\beta _{t-1}-\beta _{t-1|t-1})^{2}+Q   \\
&=F^{2}E_{t-1}P_{t-1|t-1}+Q 
\end{align}
$$
Since $E_{t-1}[\eta _{t|}]=0$ , equation (8) can be derived in the following way:
$$\begin{align}
f_{t|t-1} &=E_{t-1}(\eta _{t|t-1}^{2})=E_{t-1}(y_{t}-y_{t|t-1})^{2} \\
&=E_{t-1}(\beta _{t}+e_{t}-\beta _{t|t-1})^{2} \\
&=E_{t-1}((\beta _{t}-\beta _{t|t-1})+e_{t})^{2} \\
&=E_{t-1}(\beta _{t}-\beta _{t|t-1})^{2}+E_{t-1}(e_{t})^{2}  \\
&=P_{t|t-1}+R 
\end{align}
$$
The derivation of the updating equations relies on the following famous result in Probability Theory: 

**Theorem 1** (Proof is in the appendix) *If two random variables A and B are jointly normally distributed then the conditional on B, A is normally distribution with mean* 
$$ \mu_{A|B}=\mu_A+\frac{\sigma_{AB}}{\sigma_B^2}(b-\mu_B) \tag{12}$$
*and variance*
$$ \sigma _{A|B}^{2}=\sigma _{A}^{2}-\frac{\sigma _{AB}^{2}}{\sigma _{B}^{2}} \tag{13}$$

To derive equations (10) and (11), assume that Assume that $\beta _{t}$ and $\eta _{t|t-1}=y_{t}-y_{t|t-1}$ are jointly normally distributed. By denoting $A=\beta _{t}$ and $B=y_{t}-y_{t|t-1}=\eta_{t|t-1}$ we can obtain: 
$$
\begin{align}
\mu _{A}&=\beta _{t|t-1} \\
\mu _{B}&=E_{t-1}(\beta _{t}+e_{t}-\beta _{t|t-1}x_{t})=E_{t-1}(\beta
_{t}-\beta _{t|t-1})=0\\
\sigma _{A}^{2}&=P_{t|t-1} \\
\sigma _{B}^{2}&=f_{t|t-1} \\
\sigma _{AB} &=E_{t-1}((\beta _{t}-\beta _{t|t-1})(\eta _{t|t-1})) \\
&=E_{t-1}((\beta _{t}-\beta _{t|t-1})(y_{t}-y_{t|t-1})) \\
&=E_{t-1}((\beta _{t}-\beta _{t|t-1})(\beta _{t}+e_{t}-\beta
_{t|t-1})) \\
&=E_{t-1}((\beta _{t}-\beta _{t|t-1})((\beta _{t}-\beta
_{t|t-1})+e_{t})) \\
&=E_{t-1}(\beta _{t}-\beta _{t|t-1})^{2}+E_{t-1}(e_{t}(\beta
_{t}-\beta _{t|t-1})) \\
&=P_{t|t-1}
\end{align}
$$

Then

$$ \begin{align}
\beta _{t|t} &=\beta _{t|t-1}+\frac{P_{t|t-1}}{f_{t|t-1}}\eta _{t|t-1} \tag{14}\\
&=\beta _{t|t-1}+K_{t}\eta _{t|t-1} 
\end{align}
$$

Now due to (13):
$$\begin{align}
P_{t|t} &=P_{t|t-1}-\frac{P_{t|t-1}^{2}}{f_{t|t-1}} \tag{15}\\
&=P_{t|t-1}-\frac{P_{t|t-1}}{f_{t|t-1}}P_{t|t-1}x_{t}\\
&=P_{t|t-1}-K_{t}P_{t|t-1} 
\end{align}
$$


One important question that must be answered is what must be the initial values for $\beta_{0|0}$ and $P_{0|0}$? Kim and Nelson (1999) suggest that those can be estimated with other parameter using MLE as follows.

From (5) we have:
$$\begin{eqnarray}
\beta _{0|0} &=&\mu +F\beta _{0|0} \\
&\Downarrow & \\
\beta _{0|0} &=&\frac{\mu }{1-F}  \nonumber
\end{eqnarray}
$$


and from (6) we have:
$$\begin{eqnarray}
P_{0|0} &=&F^{2}P_{0|0}+Q \\
&\Downarrow & \\
P_{0|0} &=&\frac{Q}{1-F^{2}}
\end{eqnarray}
$$

In the next post, I will write about how to estimate this simple model using MLE and R. 


# Appendix
**Proof of Theorem 1**
If *A* and *B* are jointly normally distributed their joint density is: 
$$f(A,B) = \frac{1}{2 \pi \sigma_A \sigma_B \sqrt{1 - \rho^2}} \exp\left(-\frac{1}{2(1 - \rho^2)} \left(\frac{(a - \mu_A)^2}{\sigma_A^2} - 2 \rho \frac{(a - \mu_A)(b - \mu_B)}{\sigma_B \sigma_B}+ \frac{(b - \mu_B)^2}{\sigma_B^2} \right)\right)$$

and the marginal distribution of *B* is:
$$ g(B)=\frac{1}{\sqrt{2\pi}\sigma_B}\exp \left(-\frac{1}{2} \left(\frac{(b - \mu_B)^2}{\sigma_B^2} \right)\right) $$
To simplify notations, let
$$ v=\frac{a-\mu _{A}}{\sigma _{A}} $$
and 
$$ u=\frac{b-\mu _{B}}{\sigma _{B}} $$
Rewrite the formula of the joint density using these simplifications:
$$f(A,B) = \frac{1}{2 \pi \sigma_A \sigma_B \sqrt{1 - \rho^2}} \exp\left(-\frac{1}{2(1 - \rho^2)} \left(v^2 - 2 \rho v u + u^2 \right)\right)$$

and marginal density function of *B* as
$$ g(B)=\frac{1}{\sqrt{2\pi}\sigma_B}\exp \left(-\frac{1}{2} \left(u^2 \right)\right) $$
Now for any two random variables *A* and *B*, conditional density $w(A|B)=\dfrac{f(A,B)}{g(B)}$ .  Using the two formulas above, the conditional density formula can be rewritten as:

$$\begin{align}
w(A|B) &=\frac{\dfrac{1}{2\pi \sigma _{B}\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-
\dfrac{1}{2(1-\rho ^{2})}[u^{2}-2\rho uv+v^{2}]}}{\dfrac{1}{\sqrt{2\pi }
\sigma _{B}}e^{-\dfrac{1}{2}u^{2}}} \\
w(A|B) &=\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{
[u^{2}-2\rho uv+v^{2}]}{2(1-\rho ^{2})}+\dfrac{u^{2}}{2}} \\
w(A|B) &=\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{
[u^{2}-2\rho uv+v^{2}]}{2(1-\rho ^{2})}+\dfrac{u^{2}-u^{2}\rho ^{2}}{%
2(1-\rho ^{2})}} \\
w(A|B) &=\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{%
[u^{2}-2\rho uv+v^{2}]-u^{2}+u^{2}\rho ^{2}}{2(1-\rho ^{2})}} \\
w(A|B) &=\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{%
[v^{2}-2\rho uv+u^{2}\rho ^{2}]}{2(1-\rho ^{2})}} \\
w(A|B) &=\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{%
[v-\rho u]^{2}}{2(1-\rho ^{2})}} \\
w(A|B) &=\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{1}{2%
}\left[ \dfrac{v-\rho u}{\sqrt{(1-\rho ^{2})}}\right] ^{2}}
\end{align}$$

Now express the last equation in terms of original variables *A* and *B*:
$$\begin{eqnarray*}
w(A|B) &=&\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{1}{2%
}\left[ \dfrac{\frac{a-\mu _{A}}{\sigma _{A}}-\rho \frac{b-\mu _{B}}{\sigma
_{B}}}{\sqrt{(1-\rho ^{2})}}\right] ^{2}} \\
w(A|B) &=&\dfrac{1}{\sqrt{2\pi }\sigma _{A}\sqrt{1-\rho ^{2}}}e^{-\dfrac{1}{2%
}\left[ \dfrac{a-\{\mu _{A}+\rho \dfrac{\sigma _{A}}{\sigma _{B}}(b-\mu
_{B})\}}{\sigma _{A}\sqrt{(1-\rho ^{2})}}\right] ^{2}}
\end{eqnarray}$$

From the above, it must be clear that conditional on *B*, *A* is normally distributed with mean 
$$ \begin{eqnarray}
\mu _{A|B} &=&\mu _{A}+\rho \frac{\sigma _{A}}{\sigma _{B}}(b-\mu _{B})\\
&\Downarrow &  \nonumber \\
\mu _{A|B} &=&\mu _{A}+\frac{\sigma _{AB}}{\sigma _{B}\sigma _{A}}\frac{%
\sigma _{A}}{\sigma _{B}}(b-\mu _{B})  \nonumber \\
&\Downarrow &  \nonumber \\
\mu _{A|B} &=&\mu _{A}+\frac{\sigma _{AB}}{\sigma _{B}^{2}}(b-\mu _{B}) 
\nonumber
\end{eqnarray} $$
and variance
$$ \begin{eqnarray}
\sigma _{A|B}^{2} &=&\sigma _{A}^{2}(1-\rho ^{2}) \\
&\Downarrow &  \nonumber \\
\sigma _{A|B}^{2} &=&\sigma _{A}^{2}-\rho ^{2}\sigma _{A}^{2}  \nonumber \\
&\Downarrow &  \nonumber \\
\sigma _{A|B}^{2} &=&\sigma _{A}^{2}-\frac{\sigma _{AB}^{2}}{\sigma
_{B}^{2}\sigma _{A}^{2}}\sigma _{A}^{2}  \nonumber \\
&\Downarrow &  \nonumber \\
\sigma _{A|B}^{2} &=&\sigma _{A}^{2}-\frac{\sigma _{AB}^{2}}{\sigma _{B}^{2}}
\nonumber
\end{eqnarray}$$

*Q.E.D*

# References

Kim, C. J., & Nelson, C. R. (1999). State-space models with regime switching: classical and Gibbs-sampling approaches with applications. _MIT Press Books_, _1_.
