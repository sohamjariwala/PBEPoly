# Population balance-based rheological constitutive model
This program solves a population balanced based constitutive model for
describing the shear rheology of aggregating colloidal suspensions. The
equations are solved using method of moments with interpolative closure (MOMIC).

## Model equations
The model solves the following equations for simple shear flows and determines
the time-dependent stress and structural response.

* Population balance equation
```math
\frac{d\mu_k}{dt} = Br(\mu_k) + SA(\mu_k,\dot\gamma) + SB(\mu_k,\dot\gamma) 
```

* Viscoelastic stress equation
```math
\frac{d\sigma}{dt} = \frac{1}{\tau(\phi_a)}\left(-\sigma_{eff} +
\sigma_y(\phi_a)\mathrm{sgn}(\sigma) + \mu_s\eta(\phi_a)\dot\gamma\right)
```
* Back stress equation
```math
\frac{dA}{dt} = \dot\gamma_p(1-qA\mathrm{sgn}(\dot\gamma_p))
```
