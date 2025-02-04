# ReNLiD (Reheating Non-linear Dynamics)

Lattice simulation of interacting scalar and gauge fields in expanding space-time. Still under construction. \
Current implementation describessimulates the evolution of a scalar field $\varphi$ described by the following Lagrangian [1]:\
$$\mathcal{S} = \frac{1}{\lambda} \int d^{3}\xi d\tau \left[ \frac{1}{2} \left( \dot{\varphi} - \frac{\dot{a}}{a} \varphi \right)^{2} - \frac{\left( \nabla \varphi \right)^{2}}{2} - \frac{\varphi^{4}}{4} \right]$$


[1] [*Classical decay of inflaton*](https://arxiv.org/abs/hep-ph/9603378), S. Yu. Khlebnikov and I. I. Tkachev, Phys.Rev.Lett. 77 (1996) 219-222

## Installation

Installation is done using the `Makefile` by simply running the command `make`.