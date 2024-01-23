# Geometric and algorithmic issues in control and learning

## [PDE-AI projet kickoff meeting](https://pde-ai.math.cnrs.fr/events/kickoff_jan_24) - Jussieu, January 23-24 2024

## Nice team

[Jean-Baptiste Caillau](https://caillau.perso.math.cnrs.fr), UniCA (PR) - *Optimal control*

[Thibaud Kloczko](https://github.com/tkloczko), Inria (IR) - *Scientitic computing, dev*

[Ludovic Rifford](https://math.univ-cotedazur.fr/~rifford), UniCA (PR) - *Control, SR geometry, transport*

[Samuel Vaiter](https://samuelvaiter.com), CNRS (CR) - *Optimisation, machine learning*

## Theme 1 – The dynamics of Neural Networks training
- *Objective 1.2 – Control and machine learning: there and back again* (participants: Strasbourg, Nice) 
- **Task: Control for neural network analysis** (1 postdoc)
- *Objective 1.3 – Scalable solvers and softwares* (participants: U. Paris Cité, Nice)
- **Task: Automatic differentiation and control** (1 postdoc + 1 engineer)

## ResNets as discretised linear control systems

After Agrachev (SISSA), Sarychev (Florence), Scagliotti (TUM) *et al*:

- Agrachev, A. A.; Sarychev, A. V. [Control in the Spaces of Ensembles of Points](https://epubs.siam.org/doi/10.1137/19M1273049). *SIAM J. Control Optim.* **58** (2020), no. 3, 1579–1596
- Agrachev, A. A.; Sarychev, A. V. [Control on the Manifolds of Mappings with a View to the Deep Learning](https://link.springer.com/article/10.1007/s10883-021-09561-2). *J. Dyn. Control Syst.* **28** (2021), 989–1008
- A Scagliotti. [Deep Learning approximation of diffeomorphisms via linear-control systems](https://www.aimsciences.org/article/doi/10.3934/mcrf.2022036). *MCRF* **13** (2023), no. 3, 1226-1257

*ResNets* are compositions of nonlinear mappings

$$ \Phi = \Phi_M \circ \cdots \circ \Phi_1 $$

where $M$ is the depth of the neural network and where each building block is of the form (additional term *wrt.* non-residual networks)

$$ \Phi_l(x) = x + \sigma(W_lx+b_l). $$

This composition can be interpretated as the explicit Euler discretisation of the *Neural ODE*

$$ \dot{x}(t) = \sigma(W(t)x(t)+b(t)), $$

where $W$ and $b$ are now functions of time (continuum of layers), *controls*. Point of view developed, *e.g.* in Tabuada & Gharesifard:

- Tabuada, P.; Gharesifard, B. [Universal Approximation Power of Deep Neural Networks via Nonlinear Control Theory](https://arxiv.org/abs/2007.06007). arXiv:007.06007 (2020)

See also *constructive* approach for Lipschitz (RELU-like) activation function in Zuazua *et al*:

- Ruiz-Balet, D.; Zuazua, E. [Neural ode control for classification, approximation and transport](https://epubs.siam.org/doi/10.1137/21M1411433). *SIAM Review* **65** (2022), no. 3, 735-773

**Alternative point of view**: nonlinear in the data $x$ but *linear* in the parameters:

$$ \Phi_l(x) = x + G(x)u_l. $$

Composition now interpretated as the discretisation of 

$$ \dot{x}(t) = G(x(t))u(t) = \sum_{i=1}^m u_i(t)F_i(x(t)) $$

where the *smooth* vector fields $F_1,\dots,F_m$ are the columns of the nonlinear function $G$.

Ability to learn data (finite or *continuum*): controllability properties of the control system for *ensembles*.

## Ensembles

**Definition.** For $\Theta$ compact subset of $\mathbf{R}^n$ (set of possibly infinite indices of the data), define an ensemble as a continuous *injective* map from $\Theta$ to $\mathbf{R}^n$. Denote $\mathscr{E}_\Theta$ the set of ensembles.

**Example.** For $|\Theta| = N < \infty$ finite, ensemble = open subset of $(\mathbf{R}^n)^N$ of two by two distincts vectors: $(x_1,\dots,x_N) \in (\mathbf{R}^n)^{(N)}$ (*i.e.* $x_j=x_{j'} \implies j=j'$).

## Exact controllability

For an admissible control $u$ in $L^2([0,1],\mathbf{R}^m)$ (+ growth conditions on vector fields), define the time $1$ flow $\Phi_u$ of the controlled system

$$ \dot{x}(t) = \sum_{i=1}^m u_i(t)F_i(x(t)) \tag{1} $$

mapping an initial condition $x_0$ to $\Phi_u(x_0) := x(1,u,x_0)$.

**Definition.** The control system (1) is said to be controllable on $\mathscr{E}_\Theta$ if, for any ensembles $\gamma_0$, $\gamma_f$, there exists an admissible control $u$ *s.t.*

$$ \Phi_u \circ \gamma_0 = \gamma_f. $$

**Example.** For $|\Theta| = N$ finite and any ensembles $(x^0_1,\dots,x^0_N)$, $(x^f_1,\dots,x^0_N)$  in $(\mathbf{R}^n)^{(N)}$, controllability means that there exists an admissible control $u$ *s.t.*

$$ \Phi_u(x^0_j) = x^f_j,\quad j=1,\dots,N. $$

**Theorem.** (Agrachev-Sarychev'2020) For $k \geq 2$ controls, any $N \geq 1$ and $k$ sufficiently large, the set of vector fields $F_1,\dots,F_m$ *s.t.* controllability holds on $\mathscr{E}_\Theta$ with $|\Theta| = N$, is residual in $\mathscr{C}^k(\mathbf{R}^n)$.

**Remark.** Typical control-geometric proof: for *generic* vector fields $F_1,\dots,F_m$, the family of folds $F^{(N)}_1,\dots,F^{(N)}_m$ is bracket generating on $(\mathbf{R}^n)^{(N)}$:

$$ \text{Lie}_{x^{(N)}} \lbrace F^{(N)}_1,\dots,F^{(N)}_m \rbrace = (\mathbf{R}^n)^N,\quad x^{(N)} \in (\mathbf{R}^n)^{(N)}, $$

where the *fold* of vector field $F$ is $F^{(N)}(x_1,\dots,x_N) := (F(x_1),\dots,F(x_N)$. Controllability then ensured by Chow-Rashewsky.

## Approximate reachability

**Definition.** Let $\gamma_0$ and $\gamma_f$ be two ensembles in $\mathscr{E}_\Theta$. Then $\gamma_f$ is said to be $\mathscr{C}^0$-*approximately reachable* from $\gamma_0$ by control system (1) if, for any $\varepsilon > 0$, there exists an admissible control $u$ *s.t.*

$$ \sup_{\theta \in \Theta} |\Phi_u(\gamma_0(\theta))-\gamma_f(\theta)| \leq \varepsilon. $$

**Remark.** As a flow, any such $\Phi_u$ is diffeotopic to identity since

$$ x(0,u,\cdot) = \text{Id},\quad x(1,u,\cdot) = \Phi_u. $$

So if $\gamma_f$ is reachable from $\gamma_0$, the two ensembles must be diffeotopic.

**Definition.** The family $F_1,\dots,F_m$ satisfies the *(Lie algebra) strong approximation property* if there exists $k \geq 1$ *s.t.*, for any $\mathscr{C}^k$ vector field $X$ and for any compact $K \subset \mathbf{R}^n$, there is $\delta > 0$ *s.t.*

$$ \inf \lbrace \max_{x \in K} |X(x)-Y(x)| \text{ with } Y \in \text{Lie} \lbrace  F_1,\dots,F_m \rbrace,\ ||Y||_{1,K} \leq \delta \rbrace = 0. $$

**Remark.** The strong approximation property implies that, for every $N \geq 1$, the folds of $F_1,\dots,F_m$ are bracket generating on $(\mathbf{R}^n)^{(N)}$. So (exact) controllability holds for finite sets $\Theta$.

**Example.** The family below satisfies the strong approximation property ($n \geq 2$):

$$ F_i(x) = \frac{\partial}{\partial x_i},\quad 
  G_i(x) = e^{-|x|^2} \frac{\partial}{\partial x_i},\quad i=1,\dots,n. $$

**Theorem.** (Agrachev-Sarychev'2021) Let $\gamma_0$ and $\gamma_f$ be two *diffeotopic* ensembles in $\mathscr{E}_\Theta$. Under the strong approximation property, $\gamma_f$ is  $\mathscr{C}^0$-*approximately reachable* from $\gamma_0$ by control system (1). 

**Remark.** The proof relies on the possibility of approximating uniformly on a given compact $K$ a diffeomorphism $\Psi$ (diffeotopic to identity) by a flow $\Phi_u$ for some admissible control $u$:

$$ \sup_{x \in K} |\Phi_u(x)-\Psi(x)| \leq \varepsilon. $$

But when $\Psi$ is *not* a flow, controls cannot remain bounded when $\varepsilon \to 0$ (for instance when one triangulates $K$ with a finite number $N$ of points and uses ensemble controllability: overfitting when $N \to \infty$).

## Going further

Use optimal control to penalise / regularise according to

$$ \frac{1}{N}\sum_{j=1}^N \text{loss}( \Phi_u(x^0_j)-\Psi(x^0_j) ) + \frac{\alpha}{2} ||u||_2^2 \to \min \tag{2} $$

where $(x^0_1,\dots,x^0_N) \in (\mathbf{R}^n)^{(N)}$ is the training set of values in $K$ for the diffeomorphism $\Psi$.

- $\Gamma$-convergence results (Scagliotti'2022) when $N \to \infty$ towards

$$ \int_K \text{loss}( \Phi_u(x)-\Psi(x) )\ \mathrm{d}\mu + \frac{\alpha}{2} ||u||_2^2 $$
 
- numerical methods: gradient flow on (2), indirect method (Pontrjagin maximum principle)
- need for direct / indirect optimal control solvers / AD (automatic differentiation): [control-toolbox](https://control-toolbox.org/docs/optimalcontrol/dev/juliacon2023.html)
