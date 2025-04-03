# Code and Erratum for *On Hecke and asymptotic categories for complex reflection groups*

I collected a bit of code relevant for the paper *On Hecke and asymptotic categories for complex reflection groups*
<a href="https://arxiv.org/abs/2409.01005">https://arxiv.org/abs/2409.01005</a> on this page.

We have **Python** code below that can be run online, e.g. in <a href="https://colab.research.google.com/">Colab</a>. The code is based on a **Mathematica** file that can be downloaded on this page.
We also have a **.sage** file or **.m** file, respectively, that can be downloaded from this site (or copied from below) and you can run it with SageMath or Magma.
All files will run in the respective online calculators, see either here <a href="https://sagecell.sagemath.org/">SageMath</a> or 
<a href="http://magma.maths.usyd.edu.au/calc/">Magma</a>.

An Erratum for the paper *On Hecke and asymptotic categories for complex reflection groups* can be found at the bottom of the page.

# Contact

If you find any errors in the paper *On Hecke and asymptotic categories for complex reflection groups* **please email me**:

[dtubbenhauer@gmail.com](mailto:dtubbenhauer@gmail.com?subject=[GitHub]%web-reps)

Same goes for any errors related to this page.


# Background

The Chebyshev polynomial of the second kind $U_e=U_e(X)$ are are certain well-studied and important polynomials, see <a href="https://en.wikipedia.org/wiki/Chebyshev_polynomials">https://en.wikipedia.org/wiki/Chebyshev_polynomials</a>. 
In representation theory one renormalizes them so that they satisfy the recursion

$$U_0=1,U_1=X,\text{ and otherwise }U_{e}=XU_{e-1}+U_{e-2}.$$

Their roots, the numbers $2\cos(e\pi/n)$, are everywhere in mathematics.

The crucial observation is that this recursion matches the fusion (=tensor product) rule of $SL_2(\mathbb{C})$.
An important construction due to Koornwinder for $SL_3(\mathbb{C})$, and Eier and Lidl in general, are certain polynomials, 
the higher Chebyshev polynomials, that are attached to the fusion rules of, say, $SL_N(\mathbb{C})$. All the code below is about computing and illustrating the 
roots of the higher Chebyshev polynomials.

# The Python code: computing higher Chebyshev polynomials

The following code computes the higher Chebyshev polynomials for $N=4$:

```
from functools import lru_cache
from sympy import symbols, expand

X1, X2, X3, X4 = symbols('X1 X2 X3 X4')

# Memoization decorator to cache results
@lru_cache(maxsize=None)
def a(m1, m2, m3, m4):
    # Base cases
    if (m1, m2, m3, m4) == (0, 0, 0, 0):
        return 1
    if (m1, m2, m3, m4) == (1, 0, 0, 0):
        return X1
    if (m1, m2, m3, m4) == (0, 1, 0, 0):
        return X2
    if (m1, m2, m3, m4) == (0, 0, 1, 0):
        return X3  # Note: Mathematica has this twice with X3 and X4, assuming X3 is correct
    # Note: If you meant X4 for (0,0,1,0), replace X3 with X4 above

    # If any index is negative, return 0
    if m1 < 0 or m2 < 0 or m3 < 0 or m4 < 0:
        return 0

    # Find maximum of the indices
    max_val = max(m1, m2, m3, m4)

    # Case 1: m1 is the maximum
    if m1 == max_val:
        return (X1 * a(m1 - 1, m2, m3, m4) 
                - a(m1 - 2, m2 + 1, m3, m4) 
                - a(m1 - 1, m2 - 1, m3 + 1, m4) 
                - a(m1 - 1, m2, m3 - 1, m4 + 1) 
                - a(m1 - 1, m2, m3, m4 - 1))

    # Case 2: m4 is the maximum
    elif m4 == max_val:
        return (X4 * a(m1, m2, m3, m4 - 1) 
                - a(m1, m2, m3 + 1, m4 - 2) 
                - a(m1, m2 + 1, m3 - 1, m4 - 1) 
                - a(m1 + 1, m2 - 1, m3, m4 - 1) 
                - a(m1 - 1, m2, m3, m4 - 1))

    # Case 3: m2 is the maximum
    elif m2 == max_val:
        return (X2 * a(m1, m2 - 1, m3, m4) 
                - a(m1 + 1, m2 - 2, m3 + 1, m4) 
                - a(m1 - 1, m2 - 1, m3 + 1, m4) 
                - a(m1 + 1, m2 - 1, m3 - 1, m4 + 1) 
                - a(m1 + 1, m2 - 1, m3, m4 - 1) 
                - a(m1 - 1, m2, m3 - 1, m4 + 1) 
                - a(m1, m2 - 2, m3, m4 + 1) 
                - a(m1 - 1, m2, m3, m4 - 1) 
                - a(m1, m2 - 2, m3 + 1, m4 - 1) 
                - a(m1, m2 - 1, m3 - 1, m4))

    # Case 4: m3 is the maximum
    else:  # m3 == max_val
        return (X3 * a(m1, m2, m3 - 1, m4) 
                - a(m1, m2 + 1, m3 - 2, m4 + 1) 
                - a(m1, m2 + 1, m3 - 1, m4 - 1) 
                - a(m1 + 1, m2 - 1, m3 - 1, m4 + 1) 
                - a(m1 - 1, m2, m3 - 1, m4 + 1) 
                - a(m1 + 1, m2 - 1, m3, m4 - 1) 
                - a(m1 + 1, m2, m3 - 2, m4) 
                - a(m1 - 1, m2, m3, m4 - 1) 
                - a(m1 - 1, m2 + 1, m3 - 2, m4) 
                - a(m1, m2 - 1, m3 - 1, m4))

a(1, 1, 1, 1)
```
The result is a polynomial in four variables.

# The SageMath code: plotting higher Chebyshev polynomials

All roots of higher Chebyshev polynomials have their first coordinate in the interior of the N-cusped hypocycloid of parametric equation

$$x(\theta)=(N-1)\cos(\theta)+\cos\big((N-1)\theta\big)\text{ and }y(\theta)=(N-1)\sin(\theta)-\sin\big((N-1)\theta\big).$$

Recall that these are plane curves generated by following a point on a circle of radius 1 that rolls within a circle of radius N. Here are the pictures:

![Rolling circles](https://github.com/dtubbenhauer/nhedral/blob/main/nhedral1.png)

To be completely explicit, let $N=4$ and $e=2$. The common roots of the Chebyshev-like polynomials are:

$$(0,-1,0),(0,1,0),(\sqrt{3}i,-2,-\sqrt{3}i),(-\sqrt{3}i,-2,\sqrt{3}i),(\sqrt{3},2,\sqrt{3}),(-\sqrt{3},2,-\sqrt{3}),(e^{2i\pi/8},0,e^{14i\pi/8}),(e^{6i\pi/8},0,e^{10i\pi/8}),(e^{10i\pi/8},0,e^{6i\pi/8}),(e^{14i\pi/8},0,e^{2i\pi/8}).$$

We have the following plot of the first coordinates:

![First coordinate](https://github.com/dtubbenhauer/nhedral/blob/main/nhedral2.png)

The following SageMath code will create tex code to plot these pictures.

```
#half sum of positive roots (sage uses the realization of the root system of type A_(N-1) in R^N, and we work in the usual N-1 dimensional subspace with sum of coordinates 0) 

def rho(N):
    return vector([(N+1-2*i)/2 for i in [1..N]])

#fundamental weights 

def fund(N,i):
    return vector([1-i/N for j in [1..i]]+[-i/N for j in [i+1..N]])

#computes the set X^+(e)

def param(N,e):
    return list(filter(lambda t: sum(t) <=e,tuples([0..e],N-1)))

#computes the set V'_e

def Vpar(N,t):
    return [N*fund(N,i)*(sum([t[j-1]*fund(N,j) for j in [1..N-1]])+rho(N)) for i in [1..N-1]]

#computes the j-th Koornwinder Z-function, here q is already the root of unity needed for the Koornwinder variety V_e

def Zfun(N,e,j,gamma):
    q=exp(2*I*pi/(N*(e+N)))
    gammar=[gamma[0]]+[gamma[i]-gamma[i-1] for i in [1..N-2]] + [-gamma[N-2]]
    z=sum([q^(sum([gammar[i] for i in t])) for t in list(Combinations([0..N-1],j))])
    return z.real().simplify()+I*z.imag().simplify()

#list all points in the Koornwinder variety V_e

def Koornwinder(N,e):
    return [[Zfun(N,e,j,gamma) for j in [1..N-1]] for gamma in [Vpar(N,t) for t in param(N,e)]]

#draw the nice pictures (there is a truncation up to 3 decimal places for the points)
from collections import Counter

def KoornwinderTex(N,e):
    precision=3
    koorn=Koornwinder(N,e)
    tex="\\documentclass[crop,tikz]{standalone}\n\n\\begin{document}\n"
    if N%2 == 0:
        maxi=N/2
    else:
        maxi=(N-1)/2
    for i in [1..maxi]:
        tex+="\\begin{tikzpicture}\n"
        tex+="\\def\\a{1} \\def\\b{"+str(N/gcd(N,i))+"} \\def\\scale{"+str(binomial(N,i)/(N/gcd(N,i)))+"}\n"
        tex+="\\draw[cyan!30,very thin] ({-\\b*\\scale},{-\\b*\\scale}) grid ({\\b*\\scale},{\\b*\\scale});\n"
        tex+="\\draw[->] ({-\\b*\\scale},0) -- ({\\b*\\scale},0);\n"
        tex+="\\draw[->] (0,{-\\b*\\scale}) -- (0,{\\b*\\scale});\n"
        tex+="\\draw (0,0) circle ({\\b*\\scale});\n"
        tex+="\\draw[line width=2pt,red] plot[samples=100,domain=0:\\a*360,smooth,variable=\\t] ({\\scale*((\\b-\\a)*cos(\\t)+\\a*cos((\\b-\\a)*\\t/\\a)},{\\scale*((\\b-\\a)*sin(\\t)-\\a*sin((\\b-\\a)*\\t/\\a)});\n\n"
        k=Counter([round(float(t[i-1].real()),precision)+I*round(float(t[i-1].imag()),precision) for t in koorn]) #counts the multiplicity of the points, precision of 3 decimal places
        for x in k:
            tex+="\\draw[blue,fill=blue] ("+str(x.real())+","+str(x.imag())+") circle (4*"+str(round(0.2*k[x],1))+"pt);\n"
        tex+="\\end{tikzpicture}\n\n"
    tex+="\\end{document}"
    return tex

#Example: print LaTeX code for N=4 and e=6

print(KoornwinderTex(4,6))
```

The notation is as in the paper *On Hecke and asymptotic categories for complex reflection groups*.

# The Magma code: spectrum of graphs

We additionally need to check that certain graphs have eigenvalues being (multi)subsets of the roots of the Chebyshev polynomials (so they are in the Koornwinder variety). The 
corresponding calculations can be found in the folders on this side, ordered by the graph names. For example, in 2A3-c, the Magma code is:

```
M1:=Matrix(CyclotomicField(16),12,12,[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 
1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 
0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 
0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0]);

M2:=Matrix(CyclotomicField(16),12,12,[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 
0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 
1, 1, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 
0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 
0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0]);

M3:=Transpose(M1);

D,V:=Diagonalization([M1,M2,M3]);

[[D[j][i][i] : j in [1..3]]: i in [1..12]]
```
One can plot the graphs by copying the adjacency matrices into, for example, <a href="https://graphonline.ru/en/">GraphOnline</a>. 

# Erratum

Empty so far.
