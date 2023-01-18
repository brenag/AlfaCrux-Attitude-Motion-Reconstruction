## AlfaCrux Attitude Motion Reconstruction

## Getting Started
To run this code, you must have MATLAB and its [Aerospace Toolbox](https://www.mathworks.com/help/aerotbx/) installed.

### Setting up the environment
  - Clone the `master` branch of `AlfaCrux-Attitude-Motion-Reconstruction` from [here](https://github.com/brenag/AlfaCrux-Attitude-Motion-Reconstruction/tree/main) into `~/tmp`.
  ```
  $ cd ~/tmp
  $ git clone https://github.com/brenag/AlfaCrux-Attitude-Motion-Reconstruction.git
  $ cd AlfaCrux-Attitude-Motion-Reconstruction
  $ git checkout main
  ```

## Necessary input data

### Two-line element set for SGP4 orbit propagation

Past GP element sets [TLE](https://celestrak.org/NORAD/documentation/gp-data-formats.php) can be requested on the [Celestrak](https://celestrak.org/NORAD/archives/request.php) website. The input file has to be on the `AlfaCrux-Attitude-Motion-Reconstruction/source/` directory. The TLE file must be on its standard format (no satellite name on Line 0). 

```
1 52160U 22033D   22219.78488558  .00007175  00000+0  31352-3 0  9999
2 52160  97.3869 300.2727 0009244 201.0523 159.0335 15.22365646 19489
```
### In-orbit Measurements

To facilitate the arrangement of measurements, the data can be placed in a spreadsheet or a .csv file and then converted into a MATLAB workspace using the following [script](source/TLE/csv_to_mat.m). The measurements must be organized in columns, where each line represents a sample for an instant of time. The order in which the columns are arranged is as shown in the following table.


## Extended Kalman filter testing scheme

1. [SGP4 orbit propagation](source/sgp4.m) ([SGP4](https://celestrak.org/NORAD/documentation/gp-data-formats.php))
2. [IGRF geomagnetic field determination](source/IGRF_orbital.m) ([IGRF](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html))
3. [Cost function minimization for initial conditions](source/differential_evolution.m) ([DE](https://link.springer.com/article/10.1023/A:1008202821328))
4. [Extended Kalman Filter](source/EKF.m) ([SiriusSat-1 EKF](https://www.sciencedirect.com/science/article/abs/pii/S0094576521003957))

There are 2 different testing schemes, where one only uses in-orbit magnetometer readings and other alternative test in case of undersampled magnetometer data, using the telemetry to simulate the magnetic field measurements.

-  [Standard testing scheme](source/Main_AlfaCrux_standard.m)
-  [Undersampled testing scheme](source/Main_AlfaCrux_draft.m)
  
  
## Input data examples

Telemetry packages and respective TLE can be used to test the application. These can all be found in the `source/test/data/` directory.

## Citing

If you use the AlfaCrux Attitude Motion Reconstruction, kindly cite the following associated publication. You can obtain updates and further information about the AlfaCrux mission and also follow the LODESTAR-UnB team recent activities on our home page ([LODESTAR](https://lodestar.aerospace.unb.br/projects)).

```
@Article{app12199764,
AUTHOR = {Borges, Renato Alves and dos Santos, Andrea Cristina and Silva, William Reis and Aguayo, Leonardo and Borges, Geovany Araújo and Karam, Marcelo Monte and de Sousa, Rogério Baptista and García, Bibiano Fernández-Arruti and Botelho, Vitor Manuel de Sousa and Fernández-Carrillo, José Manuel and Lago Agra, José Miguel and Agelet, Fernando Aguado and Borges, João Vítor Quintiliano Silvério and de Oliveira, Alexandre Crepory Abbott and de Mello, Bruno Tunes and Avelino, Yasmin da Costa Ferreira and Modesto, Vinícius Fraga and Brenag, Emanuel Couto},
TITLE = {The AlfaCrux CubeSat Mission Description and Early Results},
JOURNAL = {Applied Sciences},
VOLUME = {12},
YEAR = {2022},
NUMBER = {19},
ARTICLE-NUMBER = {9764},
URL = {https://www.mdpi.com/2076-3417/12/19/9764},
ISSN = {2076-3417},
DOI = {10.3390/app12199764}
}

```
