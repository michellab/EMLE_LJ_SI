### Fitting


$$
-f_6(x_{ii})\frac{C^{\text{XDM}}_{6,ii}}{r^{6}}-f_8(x_{ii})\frac{C^{\text{XDM}}_{8,ii}}{r^{8}}-f_{10}(x_{ii})\frac{C^{\text{XDM}}_{10,ii}}{r^{10}}=-\frac{C^{\text{LJ}}_{6,ii}}{r}
$$

where $f_n$ is the Tang-Toennies damping function of order $n$, which reads:

$$
f_n({x_{ij}})=1-\left(\sum_{k=0}^n \frac{x^k_{ij}}{k!} \right)\exp(-x_{ij})
$$

The argument of the damping function is the interatomic distance divided by the average of the Slater widths of the atoms involved, i.e.:

$$
x_{ij}=\frac{r_{ij}}{(s_i+s_j)/2}
$$

