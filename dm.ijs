load 'plot'
load 'tables/csv'

NB. Data Morgana v 0.1a
NB. (c) Stichting Biomedische Rekenkamer, Biodys BV, D. Georganas
NB.
NB. ==========  Parameters can be optained from a CSV file, e.g.
NB.
NB. ┌──────┬────┐
NB. │eta_ds│0.48│
NB. └──────┴────┘
params =: readcsv'params.csv'
NB. ========== Todo ^^


NB. ======================================== Defaults for code testing


reci=: %
delta=: reci 3.69
eta_m=: reci 3.48
eta_ds=: reci 3.48
eta_dc=: reci 3.48
eta_c =: reci 42
eta_s=: reci 28
mu=: reci 79.10
tau1 =: 0.1
tau2=: 0.2
tau3=: 0.1

NB. ---------- progression to severe etc.

fsa=:0.01
fs=: 0.1 0.5 1 1.2 2.3 4.5 7.8 27.6 * fsa NB. young -> old
fca=:0.002
fc =: 0.2 0.3 1 1.8 4.7 10.6 13.6 8.7 * fca
aha=: 0.0002
alpha=: 0.1 0.3 1 3.0 10.0 45.0 120.0 505.0 * aha
fm=: 1-fs+fc


beta=: 1.52 

sigma =: 0.53 0.57 0.58 0.62 0.70 0.52 0.50 0.47
lambda=: 8 # 0.14213
ksi=:0
DEm=:0.1
DEsc=:0.1

NB. ---------- input values placeholders

inputs =:?  ,.>8# < 100000 100  
inputs =:  ,.>8# < 100000 10   , 12 # 10



NB. ================================================== Treatment for severe and critical cases


sirseci =: monad define
'S E Im Rm Is Ic Ds Dc Rsc Is_t Ic_t Ds_t Dc_t Rsc_t'  =: y
NB. -------------------- Untreated population
S1=:  (S - (lambda + mu + ksi)*S) + DEsc * (eta_ds*Is_t) + (eta_dc*Ic_t)
E1=:  (E + lambda * S) - (delta + mu + ksi)  * E
Im1=:  (Im + fm * delta * E) - (eta_m + mu + ksi ) * Im
Rm1=:  (Rm + (eta_m*Im) + (eta_s*Ds) + (eta_c*Dc)) - (mu + ksi)* Rm   NB. Rm addto
Is1=:  (Is + (fs * delta * E)) - (eta_ds + mu + ksi + tau1)*Is	      NB. add tau
Ic1=:  (Ic + (fc * delta * E)) - (eta_dc + mu + ksi + tau1)*Ic	      NB. add tau
Ds1=:  (Ds + eta_ds*Is) - (eta_s + mu + ksi)*Ds
Dc1=:  (Dc + eta_dc*Ic) - (eta_c + mu + ksi)*Dc
Rsc1=: (Rsc + (eta_s * Ds) + (eta_c*Dc)) - (mu + ksi)*Rsc
NB. -------------------- Treated population
Is_t1 =: (Is_t + (tau1 * Is))  - (eta_ds + ksi + mu) *Is_t
Ic_t1 =: (Ic_t + (tau1 * Ic))  - (eta_dc + ksi + mu) *Ic_t
Ds_t1 =: (Ds_t + ((1-DEsc)*eta_ds*Is_t)) -(eta_s+mu+ksi)*Ds_t
Dc_t1 =: (Dc_t + ((1-DEsc)*eta_dc*Ic_t)) -(eta_c+mu+ksi + alpha)*Dc_t
Rsc_t1=: (Rsc_t + (eta_s*Ds_t) + (((DEm*alpha) + eta_c)*Dc_t)) - (mu+ksi)*Rsc_t
S1, E1, Im1, Rm1, Is1, Ic1, Ds1, Dc1, Rsc1, Is_t1, Ic_t1, Ds_t1, Dc_t1,: Rsc_t1
)

prseci=: monad define
pd 'reset'
pd 'subtitle sir''r  us'
pd 'keyfont Arial 3'
pd 'key S E Im Rm Is Ic Ds Dc Rsc Is_t Ic_t Ds_t Dc_t Rsc_t' 
pd 'keypos hoc'
pd 'ycaption count'
pd 'sub 4 2'
pd 'ylog 1'
ff=: <"2 |: sirseci  ^:(i. y) |:inputs
n=.0
while. n < 8 do.
if. n=0 do.
pd 'use'
pd 'key S E Im Rm Is Ic Ds Dc Rsc Is_t Ic_t Ds_t Dc_t Rsc_t' 
pd 'keypos tr'
else.
pd 'new'
pd 'ylog 1'
end.
pd >:> n{ff
n=. >:n
end.
pd 'show'
)

NB. ================================================== Post-exposure treatment

DEpostep=:0.3
DEsc =: 0.2
DEp =: 0.3


sirpost =: monad define
'S E Im Is Ic Ds Dc R E_t Im_t Is_t Ic_t Ds_t R_t' =: y

NB. -------------------- Untreated population
S1  =: (S - (lambda + mu + ksi) *S ) + (DEpostep * delta * E_t)
E1  =: (E + (lambda * S)) - (delta + mu + ksi + tau2) * E
Im1 =: (Im + (fm * delta * E))  - (eta_m + mu + ksi) * Im
Is1 =: (Is + (fs * delta * E)) - (eta_s + mu + ksi) * Is
Ic1 =: (Ic + (fc * delta * E)) - (eta_c + mu + ksi) * Ic
Ds1 =: (Ds + (eta_ds * Is)) - (eta_s + mu + ksi) * Ds
Dc1 =: (Dc + (eta_dc * Ic)) - (eta_c + mu + ksi + alpha) * Dc
R1  =: (R  + (eta_m * Im) + (eta_s * Ds) + (eta_c * Dc)) - (mu + ksi)*R

NB. -------------------- Treated population
E_t1  =: (E_t + (tau2*E)) - (delta + mu + ksi)*E_t
Im_t1 =: Im + ((1 - DEsc +   DEsc % fm) * (1 - DEpostep) * fm * delta * E_t)  - ( mu + ksi + eta_m % 1 - DEp)*Im_t
Is_t1 =: (Is_t + ((1 - DEsc) * (1 - DEpostep) * fs * delta * E_t)) - (eta_ds * mu + ksi) * Is_t
Ic_t1 =: (Ic_t + ((1 - DEsc) * (1 - DEpostep) * fc * delta * E_t)) - (eta_dc * mu + ksi) * Ic_t
Ds_t1 =: (Ds_t + (eta_ds*Is_t)) - (eta_s + mu + ksi)*Ds_t
Dc_t1 =: (Dc_t + (eta_dc*Ic_t)) - (eta_c + mu + ksi)*Dc_t
R_t1  =: R_t + ((eta_m % 1 - DEp) * Im_t) + (eta_s * Ds_t) + (eta_c * Dc_t) - (mu + ksi)*R_t

S1, E1, Im1, Is1, Ic1, Ds1, Dc1 , R1, E_t1, Im_t1, Is_t1, Ic_t1, Ds_t1,: R_t1
)


prpost=: monad define
pd 'reset'
pd 'subtitle sir''r  us'
pd 'keyfont Arial 3'
pd 'key S E Im Rm Is Ic Ds Dc Rsc Is_t Ic_t Ds_t Dc_t Rsc_t' 
pd 'keypos hoc'
pd 'ycaption count'
pd 'sub 4 2'
NB. pd 'ylog 1'
ff=: <"2 |: sirpost  ^:(i. y) |: inputs
n=.0
while. n < 8 do.
if. n=0 do.
pd 'use'
pd 'key S E Im Rm Is Ic Ds Dc Rsc Is_t Ic_t Ds_t Dc_t Rsc_t' 
pd 'keypos tr'
else.
pd 'new'
NB. NB. pd 'ylog 1'
end.
pd >:> n{ff
n=. >:n
end.
pd 'show'
)

NB.================================================== Pre-exposure treatment


DEprep=: 0.1
DEsc=: 0.1
DEp =: 0.1


inputs2=:10 ,.~10,.~inputs NB. extended placeholders

sipre =: monad define
'S E Im Is Ic Ds Dc R T E_t Im_t Is_t Ic_t Ds_t Dc_t R_t' =: y

NB. -------------------- Untreated population
S1  =: S - (((lambda + mu + ksi) - tau3) * S)
E1  =: (E + (lambda * S))  - (delta + mu + ksi) * E
Im1 =: Im + (fm * delta * E) - (eta_m + mu + ksi) * Im
Is1 =: Im + (fs * delta * E) - (eta_ds + mu + ksi) * Im
Ic1 =: Ic + (fc * delta * E) - (eta_dc + mu + ksi) * Ic
Ds1 =: Ds + (eta_ds * Is) - (eta_s + mu + ksi) * Ds
Dc1 =: Dc + (eta_dc * Ic) - (eta_c + mu + ksi + alpha) * Dc
R1  =: R + (eta_m * Im) + (eta_s * Ds ) + (eta_c * Dc) - (mu + ksi) * R

NB. -------------------- Treated population
T1    =: T + (tau3 * S) - (((1 - DEprep) * lambda) + mu  + ksi)  * T
E_t1  =: E_t + ((1 - DEprep) * lambda * T) - (delta + mu + ksi) * E_t
Im_t1 =: Im_t  + ((1 - DEsc  + DEsc % fm) * fm * delta * E_t)  - (mu + ksi + eta_m % 1 - DEp) * Im_t
Is_t1 =: Is_t  + ((1 - DEsc) * fs * delta * E_t)  - (eta_ds + mu + ksi) * Is_t
Ic_t1 =: Ic_t  + ((1 - DEsc) * fc * delta * E_t)  - (eta_dc + mu + ksi) * Ic_t
Ds_t1 =: Ds_t  + (eta_ds * Is_t) - (eta_s + mu + ksi) * Ds_t
Dc_t1 =: Dc_t  + (eta_dc * Ic_t) - (eta_c + mu + ksi + alpha) * Dc_t
R_t1  =: R_t  +   ( Im_t *   eta_m %  1 - DEp ) + (eta_s * Ds_t) + (eta_c * Dc_t) - (mu + ksi)* R_t
S1, E1, Im1, Is1, Ic1, Ds1, Dc1, R1, T1, E_t1, Im_t1, Is_t1, Ic_t1, Ds_t1, Dc_t1,: R_t1
)


prpre =: monad define
pd 'reset'
pd 'subtitle sir''r  us'
pd 'keyfont Arial 3'
pd 'key S E Im Rm Is Ic Ds Dc Rsc Is_t Ic_t Ds_t Dc_t Rsc_t' 
pd 'keypos hoc'
pd 'ycaption count'
pd 'sub 4 2'
NB. pd 'ylog 1'
ff=: <"2 |: sipre  ^:(i. y) |: inputs2
n=.0
while. n < 8 do.
if. n=0 do.
pd 'use'
pd 'key S E Im Rm Is Ic Ds Dc Rsc Is_t Ic_t Ds_t Dc_t Rsc_t' 
pd 'keypos tr'
else.
pd 'new'
NB. NB. pd 'ylog 1'
end.
pd >:> n{ff
n=. >:n
end.
pd 'show'
)
