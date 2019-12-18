# libSchemes

**Modified Differencing Schemes for OpenFOAM**

- **GammaM** is a differencing scheme for eddy-resolving approaches (DDES,IDDES,LES), based on Gamma scheme. In the modified GammaM scheme there are two free parameters: γ and δ. The proposed GammaM scheme correctly describes the transfer of energy from large vortices to small ones at Cs = 0.15 and is stable at the values of the free parameters  γ=0.1, δ=[0.8,0.9].
- **filteredLinearM** is a differencing scheme for eddy-resolving approaches (DDES,IDDES,LES), based on filteredLinear scheme. In the modified filteredLinearM scheme, there are two free parameters α and β. The proposed filteredLinearM scheme at the values of the free parameters of α = 1 and β=[0.85,0.9] and Smagorinsky’s constant of Cs = 0.15 correctly describes the energy spectrum and is more stable than its original implementation. 
- **Koren** is a differencing scheme (flux limiter) for sufficiently data.

The library supports both OpenFOAM branches (v6.0 and v1812 or earlier).

Full investigation of modified differencing schemes presents in paper Epikhin (2019). (See below)

If you use this library in your research, please cite next articles:
> Epikhin, A. (2019), Numerical Schemes and Hybrid Approach for the Simulation of Unsteady Turbulent Flows. Mathematical Models and Computer Simulations, Vol. 11, No. 6, pp. 1019–1031.
