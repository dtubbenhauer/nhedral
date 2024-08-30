# Code and Erratum for *On Hecke and asymptotic categories for complex reflection groups*

I collected a bit of SageMath and Magma code relevant for the paper *On Hecke and asymptotic categories for complex reflection groups*
<a href="https://arxiv.org/abs/2307.03044">https://arxiv.org/abs/2307.03044</a> on this page.

The code is in a **.sage** file or **.m** file, respectively, that can be downloaded from this site (or copied from below) and you can run it with SageMath or Magma.
All files will run in the respective online calculators, see either here <a href="https://sagecell.sagemath.org/">SageMath</a> or 
<a href="http://magma.maths.usyd.edu.au/calc/">Magma</a>.

An Erratum for the paper *On Hecke and asymptotic categories for complex reflection groups* can be found at the bottom of the page.

# Contact

If you find any errors in the paper *On the asymptotic category* **please email me**:

[dtubbenhauer@gmail.com](mailto:dtubbenhauer@gmail.com?subject=[GitHub]%web-reps)

Same goes for any errors related to this page.


# Background

The Chebyshev polynomial of the second kind $U_k=U_k(X)$ are are certain well-studied and important polynomials, see <a href="https://en.wikipedia.org/wiki/Chebyshev_polynomials">https://en.wikipedia.org/wiki/Chebyshev_polynomials</a>. 
In representation theory one renormalizes them so that they satisfy the recursion

$$U_0=1,U_1=X,\text{ and otherwise }U_{k}=XU_{k-1}+U_{k-1}.$$

Their roots, the numbers $2\cos(k\pi/n)$, are everywhere in mathematics.

The crucial observation is that this recursion matches the fusion (=tensor product) rule of $SL_2(\mathbb{C})$.
An important construction due to Koornwinder for $SL_3(\mathbb{C})$, and Eier and Lidl in general, are certain polynomials, 
the Koornwinder polynomials, that are attached to the fusion rules of, say, $SL_N(\mathbb{C})$. All the code below is about computing and illustrating the 
roots of the Koornwinder polynomials.

# The SageMath code

To do

# The Magma code

To do

# Erratum

Empty so far.
