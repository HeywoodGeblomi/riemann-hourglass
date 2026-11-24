# Phase 2: formalizing the hourglass

Heywood Geblomi
21 November 2026 (America/New_York)

This note takes the five steps in order. It uses only standard mathematics and computations on the first 1000 non-trivial zeros of the zeta function (computed with mpmath, 30 decimal digits). It does not prove the Riemann Hypothesis, and it does not introduce a new equation for the zeros. Existing one-pager files and the Grok project screenshot folder were left untouched.

## 1. The average-height equation

### Definitions

These are definitions, not new theorems. They attach ordinary objects in analytic number theory to the hourglass picture.

**Lock.** The completed function
\(\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)\)
satisfies \(\xi(s)=\xi(1-s)\). That functional equation is the lock between the two sides. It pairs each zero \(\rho\) with \(1-\rho\), so zeros are symmetric about the line \(\mathrm{Re}(s)=\frac12\).

**Height on the unentropic side (the zeros).** For a non-trivial zero \(\rho\) of \(\xi\) (equivalently, of \(\zeta\)), the height is the imaginary part \(t=\mathrm{Im}(\rho)\), taken non-negative. The height measure on this side is the counting measure of those ordinates. Its distribution function is the usual zero-counting function \(N(T)\), the number of non-trivial zeros with imaginary part in \((0,T]\).

**Height on the entropic side (the primes).** For a prime power \(p^k\), the height is \(k\log p\). The height measure on this side is the Chebyshev measure, whose distribution function in the coordinate \(x=e^u\) is
\(\psi(x)=\sum_{p^k\le x}\log p\).
This is the standard weighting in the explicit formula.

**Average height.** The average density of zero heights is the smooth term
\[
\bar n(T)=\frac{1}{2\pi}\log\frac{T}{2\pi},\qquad T>2\pi.
\]
Its running integral, including the constant from the Riemann-Siegel expansion, is
\[
N_{\mathrm{sm}}(T)=\frac{T}{2\pi}\log\frac{T}{2\pi e}+\frac78.
\]
That is the waist: the average height of the zero side. On the prime side the matching smooth term is the prime number theorem, \(\psi(x)\sim x\).

**Primes as gatekeepers.** The two height measures are not independent. The explicit formula
\[
\psi_0(x)=x-\sum_{\rho}\frac{x^{\rho}}{\rho}-\log(2\pi)-\frac12\log(1-x^{-2})
\]
(at non-prime-powers, \(\psi_0=\psi\); at a jump, \(\psi_0\) is the average of the left and right limits) is the identity that holds them together. The primes, through the von Mangoldt weights, are the only arithmetic input. If they are removed, the identity collapses. That is the sense in which they gatekeep the waist.

**The average-height equation.** With these definitions the equation at the waist is
\[
N(T)=N_{\mathrm{sm}}(T)+S(T)+R(T),
\]
where \(S(T)=\frac1\pi\arg\zeta(\tfrac12+iT)\) for \(T\) not equal to an ordinate, and \(R(T)=O(1/T)\). This is the Riemann-von Mangoldt formula. It is not new. The explicit formula is the same lock, written on the prime side.

### Numerical check, first 1000 zeros

The first 1000 ordinates were computed with `mpmath.zetazero`. Landmarks: \(t_1=14.13472514\), \(t_{100}=236.52422967\), \(t_{1000}=1419.42248095\). These zeros are known to lie on the critical line, so writing \(\rho_n=\tfrac12+it_n\) is justified for this range and is not an assumption of the full Riemann Hypothesis.

\(N(T)\) jumps at each \(t_n\). The fair comparison of a continuous approximation at a jump is the midpoint: we test \(n-\tfrac12\) against \(N_{\mathrm{sm}}(t_n)\), against the more precise \(1+\theta(t_n)/\pi\) (the Riemann-Siegel theta), and against \(N_{\mathrm{sm}}\) plus the first Stirling correction \(1/(48\pi T)\).

| Approximation | RMS of \(n-\tfrac12-\)approx | Mean residual | Max abs. residual |
|---|---:|---:|---:|
| \(N_{\mathrm{sm}}\) | 0.25917 | \(-0.00017\) | 0.734 |
| \(N_{\mathrm{sm}}+1/(48\pi T)\) | 0.25917 | \(-0.00019\) | 0.734 |
| \(1+\theta(T)/\pi\) | 0.25917 | \(-0.00019\) | 0.734 |

The three residuals are essentially the same because \(N_{\mathrm{sm}}\) already tracks theta to an RMS of \(3\times 10^{-5}\), and with the \(1/(48\pi T)\) term the difference from theta drops to \(4.7\times 10^{-9}\). The leftover of size about 0.26 is the usual fluctuation \(S(T)\), not a defect of the smooth term. Near \(T=t_{1000}\), the mean spacing \(2\pi/\log(T/2\pi)\) is about 1.16, while the last actual gap was 0.73, which is the expected size.

The explicit formula, truncated to these 1000 zeros, was compared to exact \(\psi_0(x)\) at nine values of \(x\) from 20 to 10,000. The relative error stayed below \(10^{-3}\) at every tested point (for example, \(4.7\times 10^{-4}\) at \(x=10{,}000\)). That is the expected truncation error, and it confirms that the primes and the zeros are talking through the known identity, not through a new one.

**Verdict on step 1.** Once height, average height, and the gatekeeper role of the primes are defined in the standard way, the average-height equation reduces to the Riemann-von Mangoldt formula and the explicit formula. It is not a new equation.

## 2. A Hilbert-Polya operator

**Status: open.** There is still no proven self-adjoint operator, constructed from a rule that does not already use the zeros, whose eigenvalues are the ordinates \(t_n\).

What has been proposed, in brief:

- Berry-Keating (1999) suggested that a quantization of the classical Hamiltonian \(H=xp\) should have a state count matching \(N_{\mathrm{sm}}(T)\). The operator on a simple domain is not self-adjoint in a way that produces the actual zeros. Later regularizations (cutoff models, and later work by Berry-Keating, Sierra-Townsend, and others) refine the average, not a theorem that the eigenvalues are the \(t_n\).
- Connes' adelic approach produces the zeros as an absorption spectrum (missing lines), which is a different kind of object than the eigenvalues of a single Hilbert-space Hamiltonian.
- Bender-Brody-Muller (2017) proposed an operator; the proposal has been criticized as not being a well-defined self-adjoint operator independent of the zeros, and it is not accepted as a proof.

**What was computed here.**

1. A finite Hermitian likeness. The real diagonal matrix \(\mathrm{diag}(t_1,\ldots,t_{1000})\) is Hermitian and has those eigenvalues by construction. It takes the zeros as input, so it does not explain them. Any other Hermitian matrix with the same spectrum is unitarily similar to this one. This is a likeness, not an operator found from another rule.

2. A Berry-Keating-style area count. For the standard cutoff region \(\{x>0,\,p>0,\,xp\le T,\,x\ge\ell,\,p\ge 2\pi/\ell\}\), the phase-space area in units of \(2\pi\) is independent of \(\ell\) and equals
   \[
   N_{\mathrm{BK}}(T)=\frac{T}{2\pi}\log\frac{T}{2\pi e}+1.
   \]
   This is exactly \(N_{\mathrm{sm}}(T)\) except that the constant is \(1\) rather than \(7/8\). The difference is \(1/8\), the next constant in the Riemann-Siegel expansion after the area term. On the first 1000 zeros, the RMS of \(n-\tfrac12-N_{\mathrm{BK}}(t_n)\) is 0.288, versus 0.259 for \(N_{\mathrm{sm}}\), and the mean residual is \(-0.125\), which is that \(1/8\). The area count recovers the average. It does not produce \(S(T)\) or the individual \(t_n\).

**Verdict on step 2.** Likeness of the average, not a real operator.

## 3. The 2:1 geodesic ratios and the height measure

### What the ratios actually are, in the notes

From `computational_note_d6.md` and the on-box survey files in `omega-ed6-h1plus-measure-data-main/data/`:

- They are Lyapunov exponents of the Kontsevich-Zorich cocycle, on the Prym-invariant block \(H_1^+\), along special closed Teichmuller geodesics of a recorded prototype in the locus \(\Omega E_5(6)\).
- A recorded example matrix \(M_+\) has characteristic polynomial \((x^2+14x+1)(x^2+4x+1)\). The eigenvalues are \(-7\pm 4\sqrt{3}\) and \(-2\pm\sqrt{3}\). The identity \((2+\sqrt{3})^2=7+4\sqrt{3}\) holds exactly, so the ratio of log-moduli is exactly \(2:1\). After dividing by \(\log\rho(A)\), with \(A=\begin{pmatrix}31&10\\65&21\end{pmatrix}\) and \(\rho(A)=26+15\sqrt{3}\), the normalized exponents are exactly \(\pm 2/3\) and \(\pm 1/3\).
- The notes record 20 independent closed geodesics with this exact ratio, among 8095 hyperbolic maps, and 46 matrices in all. A random walk on the same generators is nearly isotropic (about \(0.43\) and \(0.43\) after scaling to sum \(6/7\)). So even inside those notes, \(2:1\) is a property of special closed orbits, not of the average.
- The known sum on this locus is \(\lambda_1^++\lambda_2^+=6/7\) (Eskin-Kontsevich-Zorich). The related \(H^-\) exponents are plus or minus \(1/7\) (Moller / EKZ). Scaling the periodic \(2:1\) ratio by the sum \(6/7\) produces the candidates \(4/7\) and \(2/7\) for the measure average. The notes mark that as a conjecture, not a theorem, and set `promote_ready` to false.

### The test (stated before looking at the alternative residuals)

The Prym ratios inform the zeta average only if a pre-registered leading-density alternative that inserts \(6/7\), \(4/7\), \(2/7\), or \(1/7\) into the coefficient of \(T\log(T/(2\pi e))\) has RMS residual, on the first 1000 zeros, less than or equal to the control \(N_{\mathrm{sm}}\) RMS. Constant-term swaps (\(7/8\) replaced by one of those fractions) are reported as a sanity check only. A shift of the constant by \(O(0.1)\) is not evidence of a connection, whatever the RMS does.

### Result

Control RMS: \(0.259\). Every leading-density alternative was far worse. The least bad of them (using \(1/7\) itself as the leading coefficient) had RMS \(59\), more than two hundred times the control. Using \(6/7\) or \(4/7\) as a leading coefficient produced residuals in the thousands. The constant-term swap \(7/8\to 6/7\) changed the RMS only from \(0.25917\) to \(0.25977\), as expected for a shift of \(0.018\), and was not treated as a pass.

**Decision: the ratios stay in the calculator.** They do not inform the average-height equation of step 1. No link was forced. The two objects live in different theories: a pair of Lyapunov exponents on a genus-4 translation surface, versus the density of zeta zeros.

## 4. The calculator's lane

The Prym work, including the discriminant-12 surface picture and the recorded D=5 survey, is a specialist calculator. It is not a Riemann Hypothesis machine, and nothing in this phase is an RH proof.

**What test it runs for the step-3 connection.** It supplies the Lyapunov numbers (the exact \(2:1\) ratio on the recorded example, the claimed count of 20, the known sum \(6/7\), the plus-or-minus \(1/7\) piece, and the candidates \(4/7\) and \(2/7\)). The connection test itself is the pre-registered comparison in section 3. The calculator does not compute zeta zeros and does not get to change \(N_{\mathrm{sm}}\).

**What was run from the box.** Exact linear algebra on the recorded example \(M_+\) and \(A\), which confirms the \(2:1\) ratio. The step-3 residual test on the first 1000 zeros. The 8095-map enumeration script and the full list of 20 geodesics are not in the on-box data (only summary counts and one worked example), so that survey was not re-run.

**A standing conflict, left standing.** A later procurement report on the box (`prym-map-procurement/REPORT.md`, with `data/d5_verification.json`) builds a certified Lanneau-Nguyen D=5 prototype and proves that this recorded matrix cannot arise as the \(H_1^+\) action of an affine diffeomorphism of that surface under a real-multiplication hypothesis. An exhaustive word search there found 0 exact \(2:1\) hits (closest ratio about 1.998). The 20-orbit count is therefore a recorded formal-monodromy claim, not a certified geometric count on that surface. That conflict does not change the step-3 decision: even if all 20 orbits are genuine, they still do not inform \(N(T)\).

## 5. One-line claim that can fail

The claim was written to the file `out/CLAIM.txt` and hashed before any further testing. The file is 472 bytes of UTF-8, including a single trailing newline.

> The average-height equation is found if, and only if, there is a minimal relation, not equivalent to the Riemann-von Mangoldt formula or the explicit formula, that is built from the heights \(t_n\) and \(k\log p\) defined in this report, that fails if the primes are removed, and that locates the first 1000 non-trivial zero ordinates to a pre-registered error; it is missed if every such candidate either reduces to one of those two known identities or fails that location test.

**SHA-256:** `d0a59770f542b002fefb880188f98de44b14f48e9d7eaa8714539899f51aa800`

Under this claim, the Riemann-von Mangoldt formula, which was verified in step 1, is a miss for "finding the equation," because it reduces to a known identity. That is the honest outcome of the formalization. No further test of the claim was run after the hash was recorded.

## Files

- `REPORT.md` (this note)
- `out/results.json` (all numerical summaries)
- `out/CLAIM.txt` and its SHA-256, as above
- `out/zeros_t.npy`, `out/plot_data.npz` (ordinates and plot arrays)
- `ONEPAGER_phase2.pdf`, `ONEPAGER_phase2.png`
- `src/compute_phase2.py`, `src/render_onepager.py`

Existing files under `/workspace/riemann_onepager/` outside `phase2/`, and `/workspace/grok-project-f4d1a687/`, were not modified.
