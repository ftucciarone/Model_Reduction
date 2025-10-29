## Proper orthogonal decomposition (POD)
The proper orthogonal decomposition is introduced here in the general context of approximation theory following Chatterjee (2000), Rivlin (1969) and Cordier and Bergmann (2003).
Obtaining a compact representation of data may be pursued with the multi-variate statistical method known as \emph{proper orthogonal decomposition} (POD). The target of the procedure is to reduce the number of intercorrelated variables to a smaller set of uncorrelated variables while retaining as much of the variation in the initial variables, that is finding a representing subspace of fixed dimension which is optimal in the sense that the error in the projection onto this subspace is minimized. This serves the twofold cause of order reduction and feature extraction of the so called *coherent structures*. Let 
$$\left\lbrace \boldsymbol{d}\left(\boldsymbol{x},t\right), \: \boldsymbol{x}\in\Omega, t\in\mathbb{R}^{+} \right\rbrace$$ 
be a set of observations of a random process over a spatial domain $\Omega$. A coherent structure, as defined by Lumley (1970), is a deterministic function $\boldsymbol{\phi}$ which is best correlated, on average, with the realizations of $\boldsymbol{d}$. In other words, the functions $\boldsymbol{\phi}$ are those functions that possess the largest mean-square projection on the observations $\boldsymbol{d}$, that is $\lvert \langle \boldsymbol{d}, \boldsymbol{\phi} \rangle \rvert^{2}$. The interest on the functions $\boldsymbol{\phi}$ is in their spatial structures, so the amplitude of these functions should not be of impact on the choice, hence they are chosen to be normalised as $\lVert \boldsymbol{\phi} \rVert^{2}=1$ and the projection itself must be normalised by the norm of the the function. One can define a subspace $S$ spanned by a set of coherent structures $`\boldsymbol{\phi}_{j}`$, with $j=1,\ldots,n$ and thus defining the projection of $\boldsymbol{d}$ onto $S$ as
$$` P_{_{S}}\boldsymbol{d} = \sum_{j=1}^{n} \dfrac{\langle\boldsymbol{d},\boldsymbol{\phi}_{j}\rangle}{\lVert\boldsymbol{\phi}_{j}\rVert^{2}}\boldsymbol{\phi}_{j} `$$
and thus the minimization of the mean square projection of $\boldsymbol{d}$ onto $S$ can be stated as\
$`
	\min_{\phi}\overline{\:\left\lVert
	\boldsymbol{d} - P_{_{S}}\boldsymbol{d}
	\right\rVert^{2}\:}^{_{X}},
`$
where the overbar means averaging in some sense. In particular, proper orthogonal decomposition is designed to minimize the number $n$ of basis functions needed in equation $\eqref{eq:pod_projection}$. 
%
%
### Mathematical formulation
In this section, the development of \workofcite{book:Holmes96}, \workofcite{thesis:Rowley2001} and \workofcite{article:Rowley2011} is followed, describing the POD procedure in the context of general Hilbert spaces. Let $\mathcal{H}$ be ah Hilbert space with inner product $\langle \cdot, \cdot \rangle_{_{\mathcal{H}}}$ and induced norm $\lVert \cdot \rVert_{_{\mathcal{H}}}$. The set of functions $\left\lbrace \boldsymbol{\phi}_{j}\left(\boldsymbol{x}\right) \in \mathcal{H} \:\colon\: j=1,\ldots,n \right\rbrace$ is defined as the one that maximise the $X-$averaged projection of $\boldsymbol{d}$ onto $\boldsymbol{\phi}$, that is
$$`
	\max_{\boldsymbol{\phi}\in\mathcal{H}}
	\dfrac{\overline{\:\lvert\langle\boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}}\rvert^{2}\:}^{_{X}}}{\lVert\boldsymbol{\phi}\rVert^{2}_{_{\mathcal{H}}}},
`$$
subject to the constraint $\lVert \boldsymbol{\phi}\lVert^{2}=1$, to close the problem. A functional $\mathcal{J}\left[\boldsymbol{\phi}\right]$ can be defined as 
$$`
	\mathcal{J}\left[\boldsymbol{\phi}\right] = \overline{\:\lvert\langle\boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}}\rvert^{2}\:}^{_{X}} - \lambda\left(\lVert\boldsymbol{\phi}\rVert^{2}_{_{\mathcal{H}}}-1\right),
`$$
including the constraint through a Lagrange multiplier. A Gateaux derivative is performed to set to zero the infinitesimal variations $\boldsymbol{\phi} + \epsilon\boldsymbol{\psi}\in\mathcal{H}$, with $\epsilon\in\mathbb{R}$, that means
\begin{align*}
	\dfrac{\mathrm{d}}{\mathrm{d}\epsilon}\mathcal{J}\left[\boldsymbol{\phi} + \epsilon\boldsymbol{\psi}\right]&=
	\left.\dfrac{\mathrm{d}}{\mathrm{d}\epsilon}\left[
	\overline{\: \langle\boldsymbol{d},\boldsymbol{\phi}+\epsilon\boldsymbol{\psi}\rangle_{_{\mathcal{H}}} \langle\boldsymbol{d},\boldsymbol{\phi}+\epsilon\boldsymbol{\psi}\rangle_{_{\mathcal{H}}}\:}^{_{X}}
	- \lambda\langle\boldsymbol{\phi}+\epsilon\boldsymbol{\psi},\boldsymbol{\phi}+\epsilon\boldsymbol{\psi}\rangle_{_{\mathcal{H}}}
	\right]\right|_{\epsilon=0}\\
&=2\overline{\: \langle\boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}}\langle\boldsymbol{d},\boldsymbol{\psi}\rangle_{_{\mathcal{H}}}\:}^{_{X}} - 2\lambda\langle\boldsymbol{\phi},\boldsymbol{\psi}\rangle_{_{\mathcal{H}}}=0.
\end{align*}
Assuming commutation is possible between the averaging $\overline{\:\cdot\:}^{_{X}}$ and the inner product $\langle \cdot, \cdot\rangle_{_{\mathcal{H}}}$ one has
$$`
	\left\langle \: \overline{ \langle \boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}} \boldsymbol{d}\:}^{_{X}} -   
	\lambda\boldsymbol{\phi} ,\boldsymbol{\psi}\right\rangle_{_{\mathcal{H}}}=0
`$$
corresponding to the eigenproblem 
$$`
	\mathcal{R}\boldsymbol{\phi} = \lambda \boldsymbol{\phi}
`$$
where $\mathcal{R}\boldsymbol{\phi} = \overline{ \langle \boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}} \boldsymbol{d}\:}^{_{X}}$. Considering the case where the Hilbert space $\mathcal{H}$ is $L^{2}$, a natural case in fluid mechanics as it represents functions with finite kinetic energy, the linear operator $\mathcal{R}\boldsymbol{\phi}$ is
\begin{align*}
	\mathcal{R}\boldsymbol{\phi} &= \overline{\: \langle\boldsymbol{d},\boldsymbol{\phi}\rangle_{_{L^{2}}} \boldsymbol{d} \:}^{_{X}} = \overline{\: \int_{\Omega}\boldsymbol{d}\left(\boldsymbol{\xi},t\right)\boldsymbol{\phi}\left(\boldsymbol{\xi}\right)\,\mathrm{d}\boldsymbol{\xi} \: \boldsymbol{d}\left(\boldsymbol{x},t\right) \:}^{_{X}}\\
	&= \int_{\Omega} \overline{\: \boldsymbol{d}\left(\boldsymbol{\xi},t\right)\boldsymbol{d}\left(\boldsymbol{x},t\right)\:}^{_{X}} \boldsymbol{\phi}\left(\boldsymbol{\xi}\right)\mathrm{d}\boldsymbol{\xi}.
\end{align*}
that is the $X-$averaged two points autocorrelation function of $\boldsymbol{d}$. The functions are orthogonal in the sense that 
$$`
	\int_{\Omega}\boldsymbol{\phi}_{j}\left(\boldsymbol{x}\right)\boldsymbol{\phi}_{k}\left(\boldsymbol{x}\right)\, \mathrm{d}\boldsymbol{x} = \delta_{jk},
`$$

It is crucial to notice that the temporally averaged quantity corresponds to the two points auto-correlation function of the signal. An orthogonal transformation is performed to project the data onto the subspace generated by the eigenvectors of the sample covariance matrix. This gives the optimal linear manifold approximating the data, in the sense that it minimizes the average squared distance between the original signal and its reduced linear representation. 

\subsection{Properties of the decomposition}
In the following, properties of the proper orthogonal decomposition are listed. The proves of these statements are now classical and can be found in the literature.
\begin{itemize}
\item For a given a bounded domain, Hilbert-Schmidt theory applies and states that the eigenproblem has a denumerable set of solutions satisfying 
$$`
	\sum_{j=1}^{d} \int_{\Omega} R_{ij}\left( \boldsymbol{x}, \boldsymbol{x}^{\prime} \right) \phi_{j}^{(n)} \left( \boldsymbol{x}^{\prime} \right) \mathrm{d}\boldsymbol{x}^{\prime} = \lambda^{(n)} \phi_{i}^{(n)}\left(\boldsymbol{x}\right)
`$$
where $\lambda^{(n)}$ and $\boldsymbol{\phi}^{(n)}_{i}$ represent respectively the eigenvalue and eigenfunction of order $n \geq 1$. Each eigenfunction is solution of the maximization problem $\eqref{eq:pod_maximization}$ with the additional constraint of being orthogonal to all previous eigenfunctions.
\item $\mathcal{R}$ can be shown to be self-adjoint and non negative, so that all eigenvalues are positive, real and converging, that is
$$`
	\lambda^{(1)} \geq \lambda^{(2)} \geq \lambda^{(3)} \geq \ldots \geq 0, \quad \text{with} \quad \sum_{n=1}^{\infty}\lambda^{(n)} < +\infty.
`$$
\item The set of eigenfunction $\boldsymbol{\phi}^{(n)}$ form a complete orthogonal set, meaning that almost every member of the set $\left\lbrace \boldsymbol{d}\left(\boldsymbol{x},t\right), \: \boldsymbol{x}\in\Omega, t\in\mathbb{R}^{+} \right\rbrace$ can be reconstructed as
$$`
	\boldsymbol{d}\left(\boldsymbol{x},t\right) = \sum_{n=1}^{\infty} \alpha^{(n)}\left( t \right) \boldsymbol{\phi}^{(n)}\left( \boldsymbol{x} \right)
`$$
where $\alpha^{(n)}$, projections of $\boldsymbol{d}$ onto $\boldsymbol{\phi}$, can be computed through the orthogonality of the eigenfunctions $\boldsymbol{\phi}$ as
$$`
	\alpha^{(n)}\left( t \right) = \langle \boldsymbol{d}, \boldsymbol{\phi} \rangle_{_{\mathcal{H}}} = \sum_{i=1}^{d}\int_{\Omega} u_{i}\left( \boldsymbol{x},t \right) \phi_{i}^{\dagger (n)}\left(\boldsymbol{x}\right) \, \mathrm{d}\boldsymbol{x}.
`$$
\item Mercer's theorem: the two points correlation tensor $R_{ij}$ can be written as a uniformly convergent series 
$$`
	R_{ij\left( \boldsymbol{x}, \boldsymbol{x}^{\prime} \right)} = \sum_{n=1}^{\infty} \lambda^{(n)} \phi_{i}^{(n)}\left( \boldsymbol{x} \right) \phi_{j}^{\dagger (n)}\left( \boldsymbol{x}^{\prime} \right).
`$$
\item Stemming from the diagonal representation of $R_{ij}$, the decomposition of $\boldsymbol{d}$ on the eigenfunctions $\boldsymbol{\phi}$ and their orthogonality, one has that 
$$`
	\overline{\alpha^{(n)}\alpha^{\dagger(m)}}^{_{X}} = \delta_{nm}\lambda^{(n)},
`$$
that means that the coefficients $\alpha^{(n)}$ are mutually uncorrelated and their mean square value are the eignevalues themselves.
\item From Mercer's theorem and orthonormality of $\boldsymbol{\phi}^{(n)}$ one can write
$$`
	\sum_{i=1}^{d}\int_{\Omega} R_{ij}\left( \boldsymbol{x}, \boldsymbol{x}^{\prime} \right) \, \mathrm{d}\boldsymbol{x} = 
	\sum_{n=1}^{\infty} \lambda^{(n)} = E,
`$$
where $E$ represents in the case of fluids with velocity field $\boldsymbol{d}$, the Turbulent Kinetic Energy (TKE) integrated over the domain $\Omega$.
\end{itemize}
%
\subsection{Algorithmic approach}
Starting from the the formalism introduced in section \ref{sect:data_decomposition}, the experimental data is organized into a matrix $\boldsymbol{D}_{i,k}$ as explained in equation \eqref{eq:data_structure}. The temporal correlation matrix $\boldsymbol{K}=\boldsymbol{D}^{\dagger}\boldsymbol{D}\in\mathbb{R}^{n_{t}\times n_{t}}$ is then computed as
$$`
	\boldsymbol{K}_{ij}=\int_{\Omega} \boldsymbol{d}^{\dagger}\left(\boldsymbol{x},t_{i}\right) \boldsymbol{d}\left(\boldsymbol{x},t_{j}\right)\,\mathrm{d}\boldsymbol{x}.
$$
The eigen-problem 
$$`\label{eq:eigen_POD}
	\boldsymbol{K} \boldsymbol{\Psi} = \boldsymbol{\Psi} \boldsymbol{\Sigma}
`$$
is then solved in order to define the temporal structures $\boldsymbol{\psi}\left(t\right)$. The corresponding spatial structures are then computed by projection of the data on these temporal structures, that means
$$`\label{eq:spatial_POD}
	\boldsymbol{\phi}_{i}\left(\boldsymbol{x}\right) = \dfrac{1}{T} \int_{T} \boldsymbol{d}\left(\boldsymbol{x},t\right) \boldsymbol{\psi}_{i} \left( t\right) \,\mathrm{d}t
`$$
%



## References
<a id="1">[1]</a> 
Chatterjee, A. (2000),
An introduction to the proper orthogonal decomposition.
Current Science 78.7, pp. 808–817 [doi](http://www.jstor.org/stable/24103957).


<a id="2">[2]</a> 
Cordier, L. and Bergmann (2003),
Post-Processing of experimental and numerical data
[doi](  ).


<a id="3">[3]</a> 
Rivlin, T. J. (1969). 
An introduction to the approximation of functions.
1st ed., Blaisdell, Book in Numerical Analysis and Computer Science, Blaisdell Pub. Co.
