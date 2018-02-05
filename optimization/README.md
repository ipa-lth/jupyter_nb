# Overview:

## optim_tools.py
Sammlung aller Tools:
- accurate_slove
- bisection_max
- bisection_min (out of date)

- Matrizen Erstellung (M, N, D, H, P, ...)
- simulate (State space simulator)
- get_Steuerungsnormalform (eigene Implementation nach Textbook)
- reverse_x_order (Matrizen variiern für python.control canonical form to "Steuerungsnormalform")


### 2010-Dissertation-Dilyana_Yankulova-cvxpy-SCS
* Hydraulisher Aktor in Steuerungsnormalform
* Gleichungen: cvxpy
* Solver: SCS
* Parameterabhängige Nebenbedingung (Sigma_Sum Matrix)
* Hauptfunktion: "accurate_solve": Verdoppelt Iterationen wenn "inaccurate" status
  * "accurate_solve" für SCS erweitern mit:
    * Timeout
    * Prüfung der constraints mittels eigenvalues
    
