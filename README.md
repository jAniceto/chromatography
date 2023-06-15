# chromatography
A suite of chromatography programms including breakthrough curve modelling, isotherm determination, parameters fitting, chromatographic peak fitting and peak deconvolution, and more...


## Breif description of each function

- **Models**:
  - `models/klinkenberg.m` - Klinkenberg model.
  - `models/transport_dispersive_ldf.m` - Transport-Dispersive Model (TDM) considering mass transfer resistence in the solid to be dominant and using the Linear Driving Force Model (LDF) approach. Uses the Matlab pdepe function to solve the system of PDEs.
  - `models/transport_dispersive_ldf_1c.m` - Same as `transport_dispersive_ldf.m` but only works for one component.
  - `models/transport_dispersive_ldf_2c.m` - Same as `transport_dispersive_ldf.m` but only works for two components.
  - `models/transport_linear_1c.m` - Transport Model (TM) with linear isotherm for single component.
  - `models/plug_flow.m` - Plug flow model (no adsorption, no mass transfer) for a single component.

- **Fitting functions:**
  - `fit_klinkenberg.m` - Fit Klinkenberg model to experimental data.
  - `fit_tdm` - Fit Transport-Dispersive Model (`transport_dispersive_ldf.m`) model to experimental data.

- **Other tools:**
  - `gaussPeakFit.m` - Fits a gaussian curve to the data provided. Works for any number of peaks.

## Chromatographic models

### Klinkenberg model (`models/klinkenberg.m`)

Klinkenberg provides an useful approximation to the analytical solution of the Convection-Dispersion model proposed by Anzelius for the case of a single solute, an initially clean bed, frontal loading and negligible axial dispersion. According to the Klinkenberg approximation the solute concentration respect to axial distance and time is given by:

$$\frac{C}{C_\mathrm{F}} \approx \frac{1}{2} [ 1 + \text{erf}( \sqrt{\tau} - \sqrt{\xi} + \frac{1}{8 \sqrt{\tau}} + \frac{1}{8 \sqrt{\xi}} ) ]$$

$$\tau = K (t - \frac{z}{u_i})$$

$$\xi = \frac{K H z}{u_i} (\frac{1 - \varepsilon_b}{\varepsilon_b})$$


### Transport-Dispersive Model (`models/transport_dispersive_ldf.m`)

`transport_dispersive_ldf` models a chromatographic experiement according to the Transport-Dispersive Model (TDM) considering mass transfer resistence in the solid to be dominant and using the Linear Driving Force Model (LDF) approach (Glauckauf and Coates, 1947). Uses pdepe functin to solve the system of partial differential equations. 

**Example:** Concentration profiles inside the column over time for a pulse injection of 3 components.

![LDF_profile](images/LDF_profile.png)



## Fitting functions

### `fit_klinkenberg`

Fit the Klinkenberg model (H and K parameters) to the data provided. 

Any number of diferent experiments ([t c] datasets) can be used by expanding the ´exp_Cfeed´ and ´exp_tc´ arrays to include additional data. Parameters are fitted to all data simultaneously.


### `fit_tdm`

Fit the Transport-Dispersive Model (`transport_dispersive_ldf.m`) to the chromatographic data provided. The isotherm parameters (`H` in the case of linear isotherm) and the linear driving force mass transfer coefficient (`KLDF`).

Any number of diferent experiments can be fitted simultaneously by expanding `exp_Cfeed` and `exp_tc` to include additional data. Parameters are fitted to all data simultaneously.


**Example 1:** Fitting of a chromatographic peak using Transport-Dispersive Model.

![peak-fit-LDF](images/peak-fit-LDF.png)

**Example 2:** Fitting of three breakthrough experiments using Transport-Dispersive Model.

![breakthough-fit-LDF](images/breakthough-fit-LDF.png)





## Other tools

### `gaussPeakFit`

gaussPeakFit fits a gaussian curve to the data provided. Works for any number of peaks.

Gaussian function : f(x) = a * exp( -(x-b)^2 / (2*c^2) ) where, 
a is the max height of the peak,
b is the position of the center of the peak (mean), and
c controls the width of the peak (standard deviation)

INPUTS:

`xy` : matrix in the [x y] format containing x, y data points

`peakSplit` : vector containing the x points around which the data is separated into ist respective peaks. For 2 peaks, peakSplit must contain 1 element. To fit a single peak peakSplit must be an empty vector.

`a` : vector cointainig the initial estimations for the max height of the peak (a) for each peak.

`b` : vector cointainig the initial estimations for the mean of the peak (b) for each peak.

`c` : vector cointainig the initial estimations for the standard deviation of the peak (c) for each peak.

OUTPUTS:

`Exitflag` : optimization status. If exitflag = 1, optimization converged 

Peak height (`a`),  Peak mean (`b`), and Peak width (`c`) : optimizaed parameters

`R` : correlation coefficient calculated using Matlab corrcoef function 

`AARD` : absolute average relative deviation

`A` : Area under peak calculated using Matlab trapz function

Figure containing the original and fited data


Example:

```matlab
Optimized parameters for peak 1: 
Exitflag = 1 
Peak height (a) = 0.5317 
Peak mean (b) = 4.4812 
Peak width (c) = 0.1914 

Fit quality for peak 1: 
Correlation coefficient (R) = 0.9974 , for a p-value of 0.0000 
AARD = 0.2760 

Area under peak 1 (A) = 12.9850 

----- 

Optimized parameters for peak 2: 
Exitflag = 1 
Peak height (a) = 0.3864 
Peak mean (b) = 5.4173 
Peak width (c) = 0.2660 

Fit quality for peak 2: 
Correlation coefficient (R) = 0.9946 , for a p-value of 0.0000 
AARD = 0.3364 

Area under peak 2 (A) = 7.8238
```

![gaussPeakFit_2peaks](images/gaussPeakFit_2peaks.png)

