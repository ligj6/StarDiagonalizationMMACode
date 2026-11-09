# Order-by-order star diagonalization of Rashba bands

Two Wolfram Mathematica notebooks calculate the analytical magnetic-field expansion of a star-unitary transformation and the resulting diagonal band Hamiltonian. They implement **direct star diagonalization** using the Cartesian Moyal product for noncommuting kinetic momenta.

The recursion is the one described in Sec. IV A, Eqs. (22)–(25), of the accompanying manuscript on coadjoint-orbit bosonization of Rashba bands. The algorithm and its implementation are explained below so that the notebooks can be used independently of the manuscript.

## Files and models

| Notebook | Model | Default `maxOrder` |
|---|---|---:|
| [StarDiagonalizationRashba.nb](StarDiagonalizationRashba.nb) | Linear Rashba coupling | `5` |
| [StarDiagonalizationQuadraticRashba.nb](StarDiagonalizationQuadraticRashba.nb) | Quadratic Rashba-type coupling | `6` |

Both files contain the same core functions. Their Hamiltonians and zero-field eigenvectors differ. Each notebook is self-contained and uses built-in Wolfram Language functions.

The code variables `x` and `y` represent **kinetic-momentum components**, $x=k_x$ and $y=k_y$. Write $k^2=x^2+y^2$. The two orbital Hamiltonians are

$$
H_{\mathrm{linear}}=
\frac{k^2}{2m}I_2+
\Delta\begin{pmatrix}
0&-\mathrm{i}(x-\mathrm{i}y)\\
\mathrm{i}(x+\mathrm{i}y)&0
\end{pmatrix},
$$

$$
H_{\mathrm{quadratic}}=
\frac{k^2}{2m}I_2+
\Delta\begin{pmatrix}
0&(x-\mathrm{i}y)^2\\
(x+\mathrm{i}y)^2&0
\end{pmatrix}.
$$

Here $m$ is the effective mass, $\Delta$ is the coupling coefficient of the selected model, and $I_2$ is the two-dimensional identity. The units of $\Delta$ differ between the linear and quadratic models. The notebooks assume real `x`, `y`, `m`, and `\[CapitalDelta]`, with $m>0$, $\Delta>0$, and $k^2>0$.

The supplied `U0` ordinarily diagonalizes `Hmat`. Its first column selects the lower zero-field band and its second column the upper band. Thus `Dcoeff[0]` has diagonal entries $k^2/(2m)\mp\Delta k$ in the linear model and $k^2/(2m)\mp\Delta k^2$ in the quadratic model.

## Quick start

The notebook headers record Mathematica 14.0. The computational sequence below was checked in Wolfram Kernel 14.0 through order $B^2$ for both models.

1. Open one notebook in Mathematica and use a fresh kernel. The first cell calls ``ClearAll["Global`*"]``; both notebooks use the same global function and variable names.
2. Evaluate the first input cell to define the Moyal-product helpers, recursion, verification routine, assumptions, and default expansion order.
3. Evaluate the model-definition cell containing `r2`, `Hmat`, and `U0`. Check the initialized matrices with

   ```mathematica
   {clean[dag[U0] . U0], clean[dag[U0] . Hmat . U0]}
   ```

   The first matrix should be the identity and the second should be diagonal.
4. For a short first calculation, set `maxOrder = 2` **after** evaluating the first cell. Then evaluate the calculation cell:

   ```mathematica
   out = StarDiagonalize[Hmat, U0, maxOrder, B];
   Ucoeff = out["Ucoeff"];
   Dcoeff = out["Dcoeff"];
   Utrunc = out["Useries"];
   Dtrunc = out["Dseries"];
   ```

5. Inspect individual coefficients or the reconstructed series. Evaluate each display expression in its own input cell:

   ```mathematica
   Table[Ucoeff[n] // MatrixForm, {n, 0, maxOrder}]
   Table[Dcoeff[n] // MatrixForm, {n, 0, maxOrder}]
   Dtrunc // MatrixForm
   ```

6. Run the verification described below. Increase `maxOrder` to obtain additional coefficients. Recompute `out` and `chk` after changing the model, assumptions, sign convention, or expansion order.

Additional scratch cells contain copied TeX expressions for presentation. They are not part of this computational sequence and can be skipped. Follow the cells identified by their contents above rather than the saved `In[...]` labels, which depend on the original evaluation history. Keep `MatrixForm` for displaying results, rather than assigning it to matrices used in subsequent calculations.

## Moyal-product convention

The default setting is

```mathematica
starSign = -I/2;
```

It gives $[x,y]_\star=-\mathrm{i}B$. All derivatives are taken with respect to the Cartesian momentum variables `x` and `y`, and matrix multiplication retains its order.

For matrix-valued symbols $A$ and $C$,

$$
A\star C=\sum_{j\geq0}B^j\mathcal C_j[A,C],
$$

$$
\mathcal C_j[A,C]=
\frac{(-\mathrm{i}/2)^j}{j!}
\sum_{\ell=0}^{j}(-1)^\ell\binom{j}{\ell}
\left(\partial_x^{j-\ell}\partial_y^\ell A\right)
\left(\partial_x^\ell\partial_y^{j-\ell}C\right).
$$

The function `starN[A, C, j]` calculates $\mathcal C_j[A,C]$, **without** the factor $B^j$. Changing `starSign` to `I/2` reverses the commutator convention.

For series $A(B)=\sum_p B^p A_p$ and $C(B)=\sum_q B^q C_q$, the coefficient of order $n$ is

$$
[B^n](A\star C)=
\sum_{p+q+j=n}\mathcal C_j[A_p,C_q].
$$

`coeffStar` performs this convolution using associations of coefficient matrices. `starSeries` returns the resulting coefficients from order zero through the requested order.

## Recursive diagonalization: manuscript Eqs. (22)–(25)

Let $H$ be the input Hermitian Hamiltonian and $N$ its matrix dimension. Start with an ordinary unitary matrix $U_0$ such that

$$
D_0=U_0^\dagger H U_0=
\operatorname{diag}(\epsilon_1,\ldots,\epsilon_N).
$$

The aim is to construct, through $M=\texttt{maxOrder}$,

$$
U_\star=\sum_{n=0}^{M}B^nU_n,
\qquad
D=\sum_{n=0}^{M}B^nD_n,
$$

with $D$ diagonal and

$$
\begin{aligned}
U_\star^\dagger\star U_\star
&=U_\star\star U_\star^\dagger
=I_N+\mathcal O(B^{M+1}),\\
U_\star^\dagger\star H\star U_\star
&=D+\mathcal O(B^{M+1}).
\end{aligned}
$$

Here $I_N$ is the identity matrix and $\dagger$ denotes Hermitian conjugation. The indices $a$ and $b$ below label the zero-field bands.

### 1. Form the known sources — Eq. (22)

Suppose $U_0,\ldots,U_{n-1}$ have been determined. Define the known partial series $U_{<n}=\sum_{r=0}^{n-1}B^rU_r$ and extract

$$
R_n=[B^n]\left(U_{<n}^\dagger\star U_{<n}-I_N\right),
\qquad
S_n=[B^n]\left(U_{<n}^\dagger\star H\star U_{<n}\right).
$$

These are exactly the sources obtained by setting the unknown $U_n$ to zero in the order-$n$ equations. They are coefficient matrices, not quantities multiplied by $B^n$. Since $n\geq1$, subtracting the identity does not change `R` in the implementation.

In the code, `Uknown` and `UdagKnown` store the known coefficients. `R` is computed with `coeffStar`; `left = starSeries[UdagKnown, Hser, n]` and a second `coeffStar` calculation produce `S`.

### 2. Impose unitarity and diagonality — Eq. (23)

Write the new correction as $U_n=U_0K_n$. At order $n$, its contribution enters without additional Moyal derivatives, because each such derivative pair would introduce another power of $B$. The conditions therefore reduce to

$$
K_n^\dagger+K_n+R_n=0,
$$

$$
\left(K_n^\dagger D_0+D_0K_n+S_n\right)_{ab}=0,
\qquad a\ne b.
$$

All nonlinear combinations of lower-order coefficients are already contained in $R_n$ and $S_n$. The remaining problem for the new coefficient is linear.

### 3. Fix the band phases — Eq. (24)

The constraints leave a diagonal phase freedom, one real parameter per band. The implementation fixes it with

$$
(K_n-K_n^\dagger)_{aa}=0.
$$

This sets the diagonal anti-Hermitian part of $K_n$ to zero at each positive order. Coefficient matrices from another phase convention can look different even when they describe the same transformed operator.

### 4. Solve and update — Eq. (25)

For nondegenerate $D_0$, the solution is

$$
(K_n)_{aa}=-\frac12(R_n)_{aa},
\qquad
(K_n)_{ab}=
\frac{(R_n)_{ab}\epsilon_b-(S_n)_{ab}}
{\epsilon_a-\epsilon_b},
\quad a\ne b.
$$

The code then stores $U_n=U_0K_n$ and the diagonal part of

$$
D_n=S_n+K_n^\dagger D_0+D_0K_n.
$$

With these coefficients available, it recomputes the sources for order $n+1$. No simultaneous solution for all expansion orders is needed.

## Outputs and verification

`StarDiagonalize[Hmat, U0, maxOrder, B]` returns an association:

| Key | Contents |
|---|---|
| `"Ucoeff"` | Association `n -> U_n`, including `n = 0` |
| `"Dcoeff"` | Association `n -> D_n`, including `n = 0` |
| `"Useries"` | Reconstructed polynomial $\sum_{n=0}^{M}B^nU_n$ |
| `"Dseries"` | Reconstructed polynomial $\sum_{n=0}^{M}B^nD_n$ |

The aliases `Utrunc` and `Dtrunc` used in the notebooks are the last two outputs. The fourth argument, called `par` inside the function, selects the symbol used to reconstruct these polynomials; its default is `B`.

Run

```mathematica
chk = VerifyStarDiagonalization[out, Hmat, maxOrder];
```

The verifier recomputes the star products from the returned transformation:

| Key | Coefficient matrices returned for orders `0` through `maxOrder` |
|---|---|
| `"UdagUminusI"` | $U_\star^\dagger\star U_\star-I_N$ |
| `"UUdagminusI"` | $U_\star\star U_\star^\dagger-I_N$ |
| `"OffDiagD"` | Off-diagonal part of $U_\star^\dagger\star H\star U_\star$ |
| `"DfullCoeff"` | Full recomputed coefficients of $U_\star^\dagger\star H\star U_\star$ |

The first three should contain zero matrices. `DfullCoeff` should also agree with the stored `Dcoeff` at every calculated order. This last comparison matters because `Dcoeff` is stored as a diagonal matrix by construction.

The following checks all four conditions at **every** calculated order, rather than only the highest-order matrices displayed in the notebook:

```mathematica
And @@ Flatten[
  Table[
    {
      chk["UdagUminusI"][r] === zeroMat[Length[Hmat]],
      chk["UUdagminusI"][r] === zeroMat[Length[Hmat]],
      chk["OffDiagD"][r] === zeroMat[Length[Hmat]],
      clean[chk["DfullCoeff"][r] - out["Dcoeff"][r]] ===
        zeroMat[Length[Hmat]]
    },
    {r, 0, maxOrder}
  ]
]
```

The expected result is `True`. The quadratic notebook also contains a separate check of `chk["OffDiagD"][2]`; that cell examines order two specifically.

## Adapting the calculation

The recursion is written for $N\times N$ matrices, although the supplied examples both have $N=2$. To use another model, supply a Hermitian `Hmat`, an ordinary unitary `U0` that diagonalizes it, and suitable assumptions for simplification.

The current implementation makes several specific choices:

- **Field-independent input Hamiltonian.** Both the generator and verifier set `Hser = <|0 -> clean[Hmat]|>`. A Hamiltonian with explicit powers of `B` requires coefficient associations for those terms in both routines. Changing the reconstruction symbol `par` does not implement this extension.
- **Nondegenerate zero-field bands.** Eq. (25) divides by `eps[[a]] - eps[[b]]`. At a band touching, the recursion requires a different treatment of the degenerate subspace. The supplied Rashba examples exclude $k=0$ and $\Delta=0$.
- **Real-parameter simplification.** `clean` uses `ComplexExpand` followed by `FullSimplify` with `$Assumptions`. Review both the assumptions and this helper when introducing complex parameters.
- **Formal field expansion.** Successful symbolic checks establish star unitarity and diagonalization through the requested order. Whether a finite truncation approximates a physical regime also depends on the band separation and momentum-gradient scales. Inverse band splittings in Eq. (25) can make a direct expansion nonuniform near degeneracy.

Higher orders require more momentum derivatives and symbolic simplification. Begin with a small `maxOrder`, verify the coefficients, and then increase the order for the calculation of interest.
