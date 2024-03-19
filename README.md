# PIRT - Processing InfraRed Thermography
<img src="logo/PIRT_color.png" alt="Alt text" title="PIRT - Processing InfraRed Thermography">

## Introduction
**PIRT** or Processing Infrared Thermography images is a Matlab toolbox for filtering snapshots, calculating heat transfer coefficients and estimating the error. The toolbox also features custom error messages to ease the image processing procedure. It has an object-based structure where the information introduced is labelled-based to make it as clear as possible. Two computation blocks can be distinguished; filter and heat transfer modules. They are independent and can be utilized by themselves.

## Structure
As previously stated, **PIRT** has an object-based architecture. As such, the toolbox allows for the creation of an indefinite number of objects for processing different groups of images. Since **PIRT** is a builder class that generates the different objects, it includes some attributes (or properties) and methods that will be saved within the different objects. The attributes are:
- ```Thot```: 2D or 3D (if temporal dependencies are included) matrix of the heated temperature measurements.
- ```Tcold```: 2D or 3D (if temporal dependencies are included) matrix of the adiabatic wall temperatures. Tcold and Thot must have the same dimensions.
- ```Calculate_filter```: a _boolean_ flag that controls the utilization of the filter module.
- ```filter_params```: _struct_ that saves the parameters associated with the filter module.
- ```CalculateHeatTransfer```: a  _boolean_ flag that controls the utilization of the heat transfer module.
- ```CalculateHeatTransferError```: a  _boolean_ flag that controls the error computation for the heat transfer module.
- ```CalculateHeatTransferErrorMethod```: saves the method that will be used to estimate the heat transfer error.
- ```HeatTransfer_params```: _struct_ that saves the parameters associated with the heat transfer module.
- ```result```: struct where the results will be saved

## Initialization
As stated, the information is introduced into **PIRT** with labels in a similar fashion to Matlab’s functions.

### Temperature maps introduction
To introduce the Thot and Tcold images, the following line would be written:
```matlab
pirt_object = PIRT("Thot", Thot_matrix ,"Tcold", Tcold_matrix )
```
The order in which the labels and the subsequent information are introduced does not matter, as long as the associated variable to a label is directly after it.
```matlab
pirt_object = PIRT("Thot", Thot_matrix , "Tcold", Tcold_matrix ) % Valid
pirt_object = PIRT("Tcold", Tcold_matrix , "Thot", Thot_matrix )  % Valid
pirt_object = PIRT("Tcold", "Thot", Tcold_matrix , Thot_matrix )  % Not Valid
```
### Filter module
The introduction of filter data is straightforward:
```matlab
pirt_object = PIRT("Filter", Filter_info)
```
Being ```Filter_info``` is a struct with the relevant information for the filtering operation. **PIRT** includes the possibility of computing several filtering operations one after the other without any limitation. However, the larger the number of filters applied the higher the error introduced. By introducing the filter data **PIRT** will automatically compute the filtering operation. **PIRT** includes four filters that can be applied to the images: POD, mPOD, sgolay32 and Weiner filters. Each introduced filter must be a struct with two fields: Filter. Type and Filter.Parameters. The Type field must include a string with the name of the filter. The Parameters field is another struct with the relevant control parameters for the associated filter. As an example of how filters are added:
```matlab
%% One filtering operation
Filter_info.Type = "POD";
Filter_info.Parameters.selection = "both";
Filter_info.Parameters.Nmod = 45;
pirt_object = PIRT("Filter", Filter_info)

%% Two POD filters, one after the other:
Filter_info(1).Type = "POD";
Filter_info(1).Parameters.selection = "both";
Filter_info(1).Parameters.Nmod = 45;
Filter_info(2).Type = "POD";
Filter_info(2).Parameters.selection = "Thot";
Filter_info(2).Parameters .Nmod = 60;
pirt_object = PIRT("Filter", Filter_info)
```
An optional input for the Filter module is ```Filter.Crop``` which is a 2x2 matrix with limiting indexes where the images are to be cut. The matrix must be introduced with the first row corresponding to the x-limit points (columns) and the second row as the y-limit points (rows). Explaining the mathematical operations behind each filter is out of the scope of this user guide. However, how each filter can be introduced will be explained.

Among all the parameters and options that could be included in the available filters, all the filters included in **PIRT** require the definition of ```Filter_info.Parameters.selection```, which identifies the ```'Thot'```, ```'Tcold'``` or ```'both'``` matrices will be filtered.
```matlab
Filter_info.Parameters.selection = "both";
```

#### POD filter
Modal filtering based on classic Proper Orthogonal Decomposition (POD). The parameters that can be introduced in the POD filter are:
- ```Filter_info.Parameters.Threshold```: (optional) energy threshold for the number of modes selection.
- ```Filter_info.Parameters.Criterion```: (optional) selects the selection method for the number of modes.
- ```Filter_info.Parameters.Nmod```: (optional) directly select the number of modes.

####  m-POD filter
Filter based on the multi-scale proper orthogonal decomposition (m-POD) algorithm proposed by  [Mendez et al. (2019, JFM)](https://doi.org/10.1017/jfm.2019.212). The parameters that can be introduced in the m-POD filter are:
- ```Filter_info.Parameters.Type```: the method by which certain bands of frequencies are eliminated. It can be Peak Removal, Frequency decoupling or both.
- ```Filter_info.Parameters.f_acq```: acquisition frequency of the images.
- ```Filter_info.Parameters.Threshold```: (optional) energy threshold for the number of modes selection.

If the peak removal method is selected, the following parameters are also required:
- ```Filter_info.Parameters.fpeaks```: frequencies to be removed.
- ```Filter_info.Parameters.w```: window size around the selected frequencies that will be eliminated.

If the frequency decoupling method is selected, the following parameters are required:
- ```Filter_info.Parameters.N_regions```: number of regions in which the frequency range will be divided.
- ```Filter_info.Parameters.f_min```: lower bound of the frequency range.
- ```Filter_info.Parameters.f_max```: upper bound of the frequency range.

If both methods are selected, the required inputs of both will be required

#### sgolay32 filter
Savitzky-Golay filter for 3D matrix and order 2 polynomial. The parameters that can be introduced in the ```sgolay32``` filter are:
- ```Filter_info.Parameters.FRAMELEN```: (optional) order of the fitted polynomial.
- ```Filter_info.Parameters.h```: (optional) scaling value or matrix.

####  wiener3 filter
Wiener 3-D adaptive noise-removal filtering. ```wiener3``` lowpass filters an intensity 3D images that has been degraded by constant power additive noise. ```wiener3``` uses a pixel-wise adaptive Wiener method based on statistics estimated from a local neighbourhood of each pixel. The parameters that can be introduced in the wiener3 filter are:
- ```Filter_info.Parameters.noise```: (optional) additive noise to the filter.
- ```Filter_info.Parameters.kernel```: (optional) kernel matrix.
 
### Heat transfer module
This module is in charge of computing the heat transfer coefficient, Nusselt number and Stanton number. It has some required inputs for its correct functioning.

#### Heat Balance
The energy balance performed in PIRT includes:

$$h = \frac{q''\_{j} - q''\_{r} - q''\_{k} - q''\_{\mathrm{unsteady}} + \sum\_{n=1}^{N} q''\_{\mathrm{custom},n}}{T\_\mathrm{w}-T\_\mathrm{aw}}$$

In the previous equation, there are four relevant heat fluxes. The joule effect heat flux $q''\_{j}$ can be computed as $q''\_{j}=V I/A\_{board}$. After the heat flux due to the joule effect, there are three dissipation heat fluxes that are taken into account. $q''\_{r}$ represents the heat flux due to radiation, which is computed assuming that the surroundings act as a black body with a temperature equal to the one of the freestream. Therefore, it can be computed as $q''\_{r}=\sigma \epsilon \left(T\_{w}^4 - T\_{\inf}^4\right)$. $\epsilon$ represents the emissivity of the material (known due to the coating paint) and $\sigma$ is the Boltzmann constant. $q''\_{k}$ is the tangential heat flux in the thin foil sensor. It can be computed as $q''\_{k} = kt\left(\frac{\partial^2T}{\partial^2x}+\frac{\partial^2T}{\partial^2y}\right)$, with $t$ as the thickness of the TFS. This formula will be used if the HFS is selected to be a foil, if the PCB module is selected, the formula will in turn be $q''\_{k} = t\left(\lambda_x\frac{\partial^2T}{\partial^2x}+\lambda_y\frac{\partial^2T}{\partial^2y}\right)$, being $\lambda_x$ and $\lambda_y$ the thermal conductivity of the PCB in the x and y directions (the interested reader is referred to [Torre et al. (2018, IJHMT)](https://doi.org/10.1016/j.ijheatmasstransfer.2018.06.106)). Lastly, temporal unsteady losses are also considered as $q''_{\mathrm{unsteady}} = \rho c\_{P} t \frac{\partial T}{\partial t}$

#### Inputs
This module is in charge of computing the heat transfer coefficient, Nusselt number and Stanton number. It has some required inputs for its correct functioning. The label CalculateHeatTransfer must be introduced, followed by ```'Nu'``` (in order to calculate the Nusselt number), ```'h'``` (to calculate the heat transfer coefficient) or ```'St'``` (to calculate the Stanton number. Any combination of the selection labels is also allowed.
```matlab
pirt_object = PIRT('CalculateHeatTransfer','h') %Compute h
pirt_object = PIRT('CalculateHeatTransfer','Nu') %Compute Nu
pirt_object = PIRT('CalculateHeatTransfer','St') %Compute St
pirt_object = PIRT('CalculateHeatTransfer','h','Nu','St') %Compute the three options
```
Two more inputs are required for the computation of the heat transfer. The information of the material in which the thermal images were taken, introduced as ```'HTF'``` (heat foil sensor). The free stream conditions of the flow are also required.
```matlab
pirt_object = PIRT('HFS', HFS_struc, 'Conditions', Conditions_struc)
```
HFS is a struct that must contain the following fields:
- ```HFS.s```: the thickness of the thin foil sensor.
- ```HFS.rho```: the density of the thin foil sensor material.
- ```HFS.cp```: the cp of the thin foil sensor material.
- ```HFS.epsilon```: the emissivity of the thin foil sensor material’s face that is oriented towards the infrarred camera.
- ```HFS.Area```: the area of the thin foil sensor captured by the images.
- ```HFS.k```: the thermal conductivity of the thin foil sensor material.
- ```HFS.Type```: the type of HFS that was used in the experiment, it can either be ```"Foil"``` or ```"PCB"```. It will only affect if the selected choice was ```"PCB"``` and the tangential heat is computed. Optional input, it will be selected as Foil by default.
- ```HFS.lambdax```: the thermal conductivity of the thin foil sensor material in the x direction in the PCB model. Only required if the PCB model is selected.
- ```HFS.lambday```: the thermal conductivity of the thin foil sensor material in the y direction in the PCB model. Only required if the PCB model is selected.

Instead of introducing ```HFS.Area```, if the sensor is a rectangle, ```HFS.H``` and ```HFS.W``` can be introduced as the height and width of the sensor respectively. Lastly, the free stream conditions must also be included as a struct with the following parameters:
- ```Conditions.L```: characteristic length of the problem.
- ```Conditions.V```: voltage supplied to the TFS.
- ```Conditions.I```: current supplied to the TFS.
- ```Conditions.Tamb```: vector containing first the cold free stream temperature and second the hot free stream temperature.
- ```Conditions.Uinf```: free stream velocity. Only required for computing the Stanton number.

**PIRT** also includes spatial and temporal dependencies of heat transfer. If they are to be included, the flags ’TimeDer’ and ’SpatialDer’ have to be included. Some of the previously defined inputs are only required for these computations, the error system embedded in **PIRT** will serve as a guide of what is required.

##### Custom Heat balance terms
In case there is a heat-flux term that is not taken into account by **PIRT**, the user can introduce custom terms. These are introduced in a similar manner to the ```Nu```, ```St``` or ```h``` selection. The label that must be introduced is ```'CustomQ'```. After this label, _cell_ (even if it only has 1 element) containing in each element a $q''_{\mathrm{custom},n}$ term must be introduced. Each element of the ```'CustomQ'``` cell can be:
- Scalar: A scalar value that will be added to the whole heat balance element-wise.
- 2D matrix: matrix with the same dimensions as the Thot and Tcold snapshots. If the snapshots are 3D, the CustomQ element must have the same size as the first two dimensions of the snapshots and it will be added element wise in the third dimension.
- 3D matrix with the same dimensions as the Thot and Tcold snapshots.

#### Error estimation
**PIRT** includes a module for estimating the error in heat transfer computations. The flag CalculateHeatTransferError must be introduced for its calculation. This label can be followed by ```'Moffat'``` or ```'Montecarlo'``` to specify which error estimation method is used. If this method is not specified, ```'Moffat'``` will be used as default.
```matlab
pirt_object = PIRT('CalculateHeatTransferError') %Compute with the default Moffat method
pirt_object = PIRT('CalculateHeatTransferError','Montecarlo') %Compute with the Montecarlo method
```
The errors for the TFS and conditions must be included. They are specified as relative errors in such a way that the actual error is ```error = variable·variable_error```. For the Montecarlo error estimation, the nominal value is taken as the mean and the error (following the above description) is taken as the standard deviation in a normal distribution.
```matlab
pirt_object = PIRT('Error’, Error_struc) %Introduce the error information
```
The required input for this module is a struct labelled as ```'Error'```. With the following required fields:
- ```Error.errorThot```: relative error in the Thot temperature acquisition. It can be a single value or a matrix with a relative error foe earch pixel.
- ```Error.errorTcold```: relative error in the Tcold temperature acquisition. It can be a single value or a matrix with a relative error foe each pixel.
- ```Error.errorTamb```: relative error in the ambient temperature.
- ```Error.errorV```: relative error in the supplied boltage.
- ```Error.errorI```: relative error in the supplied current.
- ```Error.errorEpsilon```: relative error in the emissivity of the TFS.
- ```Error.errorrho```: relative error in the density of the TFS material.
- ```Error.errorcp```: relative error in the cp of the TFS’s material.
- ```Error.errors```: relative error in the TFS’s thickness.
- ```Error.errorA```: relative error in the TFS’s area.
- ```Error.errorkplate```: relative error in the thermal conductivity of the TFS.
- ```Error.errorLchar```: relative error in the characteristic length.
- ```Error.errork```: relative error in the air’s thermal conductivity.
- ```Error.errorUinf```: relative error in the free stream velocity. Only required if the Stanton number is computed in the heat transfer module.
