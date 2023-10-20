The Cauchy problem for a system of equations of motion of a material point in a potential field U(x) is considered:
```math
dx(t)/dt = v,              x(0) = x_0,
dv(t)/dt = - 1/m*dU/dx,    v(0) = v_0.
```
Considered numerical methods of solution
1. The Euler method with a constant step.
2. An explicit two-step Adams scheme.
3. The Runge-Kutta method of the 4th order.

Cauchy problem for a system of equations of motion
Bringing the system to a dimensionless form

```math
U(x) = 10 ^-18*cos(x^2), m = 9.1 10^-31, x_0 = 0, v_0 = 5 10^5;
So we have system in dimensionless form:
dx/dt = -v,           x(0) = 0
dv/dt = 3x^2 + 2x,    v(0) = -2 * \sqrt(9.1) / 100
```
Graphs (x(t)) constructed using the Euler method, with increased accuracy, with different steps and Maple:

![image](https://github.com/NoPainNoGane/numerical_sol_for_eq/assets/64308897/976990c0-844e-4ab6-8c54-54193d8a2005)

Graphs (v(t)) constructed using the Euler method, with increased accuracy, with different steps and Maple:

![image](https://github.com/NoPainNoGane/numerical_sol_for_eq/assets/64308897/0e08b017-1ac5-4659-bf83-a9c682539408)

Graphs (x(t)) constructed using the Adams method with different steps and Maple.

![image](https://github.com/NoPainNoGane/numerical_sol_for_eq/assets/64308897/44494a93-2220-446d-bc19-b901228ca3c2)

Graphs (v(t)) constructed using the Adams method with different steps and Maple.

![image](https://github.com/NoPainNoGane/numerical_sol_for_eq/assets/64308897/3a996856-0bf3-4375-83d4-5991aeef4116)

Graphs (x(t)) constructed using the 4th order Runge-Kutta method with different steps and Maple.

![image](https://github.com/NoPainNoGane/numerical_sol_for_eq/assets/64308897/63543bf2-df11-43b5-9e30-52478b2d3c37)

Graphs (v(t)) constructed using the 4th order Runge-Kutta method with different steps and Maple.

![image](https://github.com/NoPainNoGane/numerical_sol_for_eq/assets/64308897/21586cf6-4b94-4bf5-89f5-6c397861b8c9)




