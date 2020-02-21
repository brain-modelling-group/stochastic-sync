The model equations are divided into two parts:
1. deterministic (embedded in odefun.m)
2. stochastic (embedded in sdefun1.m)

The model is then simulated by the main function kuramoto_simulator.m.
The function numerically solves the solutions via two schemes (deterministic Euler and stochastic Euler-Maruyama) depending whether you want a model with or without noise (choice is automatic based on the noise strength parameter). This avoids using the stochastic solver when not needed, effectively speeding up the calculations.  
