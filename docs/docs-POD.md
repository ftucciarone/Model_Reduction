## Proper orthogonal decomposition (POD)
The proper orthogonal decomposition is introduced here in the general context of approximation theory following \workofcite{article:Chatterjee2000}, \workofcite{book:Rivlin1981} and \workofcite{article:Cordier2003}.
Obtaining a compact representation of data may be pursued with the multi-variate statistical method known as \emph{proper orthogonal decomposition} (POD). The target of the procedure is to reduce the number of intercorrelated variables to a smaller set of uncorrelated variables while retaining as much of the variation in the initial variables, that is finding a representing subspace of fixed dimension which is optimal in the sense that the error in the projection onto this subspace is minimized. This serves the twofold cause of order reduction and feature extraction of the so called \emph{coherent structures}. Let $$\left\lbrace \boldsymbol{d}\left(\boldsymbol{x},t\right), \: \boldsymbol{x}\in\Omega, t\in\mathbb{R}^{+} \right\rbrace$$ be a set of observations of a random process over a spatial domain $`\Omega`$. A coherent structure, as defined by \workofcite{book:Lumley1970}, is a deterministic function $\boldsymbol{\phi}$ which is best correlated, on average, with the realizations of $\boldsymbol{d}$. In other words, the functions $\boldsymbol{\phi}$ are those functions that possess the largest mean-square projection on the observations $\boldsymbol{d}$, that is $\lvert \langle \boldsymbol{d}, \boldsymbol{\phi} \rangle \rvert^{2}$. The interest on the functions $\boldsymbol{\phi}$ is in their spatial structures, so the amplitude of these functions should not be of impact on the choice, hence they are chosen to be normalised as $\lVert \boldsymbol{\phi} \rVert^{2}=1$ and the projection itself must be normalised by the norm of the the function. One can define a subspace $S$ spanned by a set of coherent structures $\boldsymbol{\phi}_{j}$, with $j=1,\ldots,n$ and thus defining the projection of $\boldsymbol{d}$ onto $S$ as
\begin{equation}\label{eq:pod_projection}
	P_{_{S}}\boldsymbol{d} = \sum_{j=1}^{n} \dfrac{\langle\boldsymbol{d},\boldsymbol{\phi}_{j}\rangle}{\lVert\boldsymbol{\phi}_{j}\rVert^{2}}\boldsymbol{\phi}_{j}
\end{equation}
and thus the minimization of the mean square projection of $\boldsymbol{d}$ onto $S$ can be stated as 
\begin{equation}
	\min_{\phi}\overline{\:\left\lVert
	\boldsymbol{d} - P_{_{S}}\boldsymbol{d}
	\right\rVert^{2}\:}^{_{X}},
\end{equation}
where the overbar means averaging in some sense. In particular, proper orthogonal decomposition is designed to minimize the number $n$ of basis functions needed in equation $\eqref{eq:pod_projection}$. 
%
%
### Mathematical formulation
In this section, the development of \workofcite{book:Holmes96}, \workofcite{thesis:Rowley2001} and \workofcite{article:Rowley2011} is followed, describing the POD procedure in the context of general Hilbert spaces. Let $\mathcal{H}$ be ah Hilbert space with inner product $\langle \cdot, \cdot \rangle_{_{\mathcal{H}}}$ and induced norm $\lVert \cdot \rVert_{_{\mathcal{H}}}$. The set of functions $\left\lbrace \boldsymbol{\phi}_{j}\left(\boldsymbol{x}\right) \in \mathcal{H} \:\colon\: j=1,\ldots,n \right\rbrace$ is defined as the one that maximise the $X-$averaged projection of $\boldsymbol{d}$ onto $\boldsymbol{\phi}$, that is
\begin{equation}\label{eq:pod_maximization}
	\max_{\boldsymbol{\phi}\in\mathcal{H}}
	\dfrac{\overline{\:\lvert\langle\boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}}\rvert^{2}\:}^{_{X}}}{\lVert\boldsymbol{\phi}\rVert^{2}_{_{\mathcal{H}}}},
\end{equation}
subject to the constraint $\lVert \boldsymbol{\phi}\lVert^{2}=1$, to close the problem. A functional $\mathcal{J}\left[\boldsymbol{\phi}\right]$ can be defined as 
\begin{equation}
	\mathcal{J}\left[\boldsymbol{\phi}\right] = \overline{\:\lvert\langle\boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}}\rvert^{2}\:}^{_{X}} - \lambda\left(\lVert\boldsymbol{\phi}\rVert^{2}_{_{\mathcal{H}}}-1\right),
\end{equation}
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
\begin{equation}
	\left\langle \: \overline{ \langle \boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}} \boldsymbol{d}\:}^{_{X}} -   
	\lambda\boldsymbol{\phi} ,\boldsymbol{\psi}\right\rangle_{_{\mathcal{H}}}=0
\end{equation}
corresponding to the eigenproblem 
\begin{equation}
	\mathcal{R}\boldsymbol{\phi} = \lambda \boldsymbol{\phi}
\end{equation}
where $\mathcal{R}\boldsymbol{\phi} = \overline{ \langle \boldsymbol{d},\boldsymbol{\phi}\rangle_{_{\mathcal{H}}} \boldsymbol{d}\:}^{_{X}}$. Considering the case where the Hilbert space $\mathcal{H}$ is $L^{2}$, a natural case in fluid mechanics as it represents functions with finite kinetic energy, the linear operator $\mathcal{R}\boldsymbol{\phi}$ is
\begin{align*}
	\mathcal{R}\boldsymbol{\phi} &= \overline{\: \langle\boldsymbol{d},\boldsymbol{\phi}\rangle_{_{L^{2}}} \boldsymbol{d} \:}^{_{X}} = \overline{\: \int_{\Omega}\boldsymbol{d}\left(\boldsymbol{\xi},t\right)\boldsymbol{\phi}\left(\boldsymbol{\xi}\right)\,\mathrm{d}\boldsymbol{\xi} \: \boldsymbol{d}\left(\boldsymbol{x},t\right) \:}^{_{X}}\\
	&= \int_{\Omega} \overline{\: \boldsymbol{d}\left(\boldsymbol{\xi},t\right)\boldsymbol{d}\left(\boldsymbol{x},t\right)\:}^{_{X}} \boldsymbol{\phi}\left(\boldsymbol{\xi}\right)\mathrm{d}\boldsymbol{\xi}.
\end{align*}
that is the $X-$averaged two points autocorrelation function of $\boldsymbol{d}$. The functions are orthogonal in the sense that 
\begin{equation}
	\int_{\Omega}\boldsymbol{\phi}_{j}\left(\boldsymbol{x}\right)\boldsymbol{\phi}_{k}\left(\boldsymbol{x}\right)\, \mathrm{d}\boldsymbol{x} = \delta_{jk},
\end{equation}

It is crucial to notice that the temporally averaged quantity corresponds to the two points auto-correlation function of the signal. An orthogonal transformation is performed to project the data onto the subspace generated by the eigenvectors of the sample covariance matrix. This gives the optimal linear manifold approximating the data, in the sense that it minimizes the average squared distance between the original signal and its reduced linear representation. 

\subsection{Properties of the decomposition}
In the following, properties of the proper orthogonal decomposition are listed. The proves of these statements are now classical and can be found in the literature.
\begin{itemize}
\item For a given a bounded domain, Hilbert-Schmidt theory applies and states that the eigenproblem has a denumerable set of solutions satisfying 
\begin{equation}
	\sum_{j=1}^{d} \int_{\Omega} R_{ij}\left( \boldsymbol{x}, \boldsymbol{x}^{\prime} \right) \phi_{j}^{(n)} \left( \boldsymbol{x}^{\prime} \right) \mathrm{d}\boldsymbol{x}^{\prime} = \lambda^{(n)} \phi_{i}^{(n)}\left(\boldsymbol{x}\right)
\end{equation}
where $\lambda^{(n)}$ and $\boldsymbol{\phi}^{(n)}_{i}$ represent respectively the eigenvalue and eigenfunction of order $n \geq 1$. Each eigenfunction is solution of the maximization problem $\eqref{eq:pod_maximization}$ with the additional constraint of being orthogonal to all previous eigenfunctions.
\item $\mathcal{R}$ can be shown to be self-adjoint and non negative, so that all eigenvalues are positive, real and converging, that is
\begin{equation}
	\lambda^{(1)} \geq \lambda^{(2)} \geq \lambda^{(3)} \geq \ldots \geq 0, \quad \text{with} \quad \sum_{n=1}^{\infty}\lambda^{(n)} < +\infty.
\end{equation}
\item The set of eigenfunction $\boldsymbol{\phi}^{(n)}$ form a complete orthogonal set, meaning that almost every member of the set $\left\lbrace \boldsymbol{d}\left(\boldsymbol{x},t\right), \: \boldsymbol{x}\in\Omega, t\in\mathbb{R}^{+} \right\rbrace$ can be reconstructed as
\begin{equation}
	\boldsymbol{d}\left(\boldsymbol{x},t\right) = \sum_{n=1}^{\infty} \alpha^{(n)}\left( t \right) \boldsymbol{\phi}^{(n)}\left( \boldsymbol{x} \right)
\end{equation}
where $\alpha^{(n)}$, projections of $\boldsymbol{d}$ onto $\boldsymbol{\phi}$, can be computed through the orthogonality of the eigenfunctions $\boldsymbol{\phi}$ as
\begin{equation}
	\alpha^{(n)}\left( t \right) = \langle \boldsymbol{d}, \boldsymbol{\phi} \rangle_{_{\mathcal{H}}} = \sum_{i=1}^{d}\int_{\Omega} u_{i}\left( \boldsymbol{x},t \right) \phi_{i}^{\dagger (n)}\left(\boldsymbol{x}\right) \, \mathrm{d}\boldsymbol{x}.
\end{equation}
\item Mercer's theorem: the two points correlation tensor $R_{ij}$ can be written as a uniformly convergent series 
\begin{equation}
	R_{ij\left( \boldsymbol{x}, \boldsymbol{x}^{\prime} \right)} = \sum_{n=1}^{\infty} \lambda^{(n)} \phi_{i}^{(n)}\left( \boldsymbol{x} \right) \phi_{j}^{\dagger (n)}\left( \boldsymbol{x}^{\prime} \right).
\end{equation}
\item Stemming from the diagonal representation of $R_{ij}$, the decomposition of $\boldsymbol{d}$ on the eigenfunctions $\boldsymbol{\phi}$ and their orthogonality, one has that 
\begin{equation}
	\overline{\alpha^{(n)}\alpha^{\dagger(m)}}^{_{X}} = \delta_{nm}\lambda^{(n)},
\end{equation}
that means that the coefficients $\alpha^{(n)}$ are mutually uncorrelated and their mean square value are the eignevalues themselves.
\item From Mercer's theorem and orthonormality of $\boldsymbol{\phi}^{(n)}$ one can write
\begin{equation}
	\sum_{i=1}^{d}\int_{\Omega} R_{ij}\left( \boldsymbol{x}, \boldsymbol{x}^{\prime} \right) \, \mathrm{d}\boldsymbol{x} = 
	\sum_{n=1}^{\infty} \lambda^{(n)} = E,
\end{equation}
where $E$ represents in the case of fluids with velocity field $\boldsymbol{d}$, the Turbulent Kinetic Energy (TKE) integrated over the domain $\Omega$.
\end{itemize}
%
\subsection{Algorithmic approach}
Starting from the the formalism introduced in section \ref{sect:data_decomposition}, the experimental data is organized into a matrix $\boldsymbol{D}_{i,k}$ as explained in equation \eqref{eq:data_structure}. The temporal correlation matrix $\boldsymbol{K}=\boldsymbol{D}^{\dagger}\boldsymbol{D}\in\mathbb{R}^{n_{t}\times n_{t}}$ is then computed as
\begin{equation}\label{eq:POD_correlation}
	\boldsymbol{K}_{ij}=\int_{\Omega} \boldsymbol{d}^{\dagger}\left(\boldsymbol{x},t_{i}\right) \boldsymbol{d}\left(\boldsymbol{x},t_{j}\right)\,\mathrm{d}\boldsymbol{x}.
\end{equation}
The eigen-problem 
\begin{equation}\label{eq:eigen_POD}
	\boldsymbol{K} \boldsymbol{\Psi} = \boldsymbol{\Psi} \boldsymbol{\Sigma}
\end{equation}
is then solved in order to define the temporal structures $\boldsymbol{\psi}\left(t\right)$. The corresponding spatial structures are then computed by projection of the data on these temporal structures, that means
\begin{equation}\label{eq:spatial_POD}
	\boldsymbol{\phi}_{i}\left(\boldsymbol{x}\right) = \dfrac{1}{T} \int_{T} \boldsymbol{d}\left(\boldsymbol{x},t\right) \boldsymbol{\psi}_{i} \left( t\right) \,\mathrm{d}t
\end{equation}
%
\subsection{Noise ansatz}
\begin{figure}[t]
	\centering
	\includegraphics[width=\imagewidth]{/Users/Francesco/Documents/These_Tucciarone_all/These_Tucciarone/03_Noise_Models/Figures/poster_POD_alpha0.png}
	\caption{Outline of the POD noise generation procedure}\label{fig:poster_POD}
\end{figure}%
Employing POD on a set of velocity fluctuations of type $\eqref{eq:time_fluctuations}$ or $\eqref{eq:en_scale_time_fluctuations}$ produces a set of velocity modes $\lbrace \boldsymbol{\phi}_{j}\left( \boldsymbol{x} \right), \lambda_{j}, j=1,\ldots,N\rbrace$ that can be used to define the noise ansatz as
\begin{equation}\label{eq:pod_noise_ansatz}
	\boldsymbol{\sigma}\left( \boldsymbol{x} \right) \mathrm{d}\mathbf{B}_{t} = \sqrt{\tau}\sum_{k=1}^{N} \lambda^{1/2}_{k} \boldsymbol{\phi}_{k}\left( \boldsymbol{x} \right) \mathrm{d}\beta^{k}_{t}
\end{equation}
with associated variance tensor computed as
\begin{equation}\label{eq:pod_variance_ansatz}
	\boldsymbol{a}\left( \boldsymbol{x} \right) = \tau \sum_{k=1}^{N} \lambda_{k}\boldsymbol{\phi}_{k}\left( \boldsymbol{x} \right) \boldsymbol{\phi}^{_{\mathrm{T}}}_{k}\left( \boldsymbol{x} \right).
\end{equation}
If a non centred noise of type $\eqref{eq:NonCentred_LUnoise}$ is considered favourable, the time average $\overline{\boldsymbol{u}}^{\, t}$ that was removed from the initial data can be re-inserted in the simulation through the Girsanov correction, defining thus the noise as
\begin{equation}\label{eq:pod_girsenov_noise_ansatz}
	\boldsymbol{\sigma}\left( \boldsymbol{x} \right) \mathrm{d}\mathbf{B}_{t} = - \overline{\boldsymbol{u}}^{\, t}\left( \boldsymbol{x} \right)\mathrm{d}t + \sqrt{\tau}\sum_{k=0}^{N} \lambda^{1/2}_{k}\boldsymbol{\phi}_{k}\left( \boldsymbol{x} \right) \mathrm{d}\beta^{k}_{t},
\end{equation}
where $\boldsymbol{\sigma}_t\mathbf{Y}_{t}$ in \eqref{eq:NonCentred_LUnoise} is defined as $\overline{\boldsymbol{u}}^{\, t}$.
