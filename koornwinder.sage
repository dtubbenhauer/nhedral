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
