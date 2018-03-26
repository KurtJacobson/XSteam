*XSteam was developed by Magnus Holmgren in 2006 and was published at www.x-eng.com.
Unfortunately that site has been unavailable for several years now, making it
very difficult to obtain the various versions of XSteam. Having used XSteam
extensively in the past, I have created this repo to ensure that it remains available.*

# XSteam

XSteam is a implementation of the IAPWS IF97 standard formulation. It
provides accurate data for water and steam and mixtures of water and steam
properties from 0 - 1000 bar and from 0 - 2000°C.
It is available for MATLAB, Excel, OpenOffice, and as a DLL.

Provided thermodynamic properties are:

  * Temperature
  * Pressure
  * Enthalpy
  * Specific volume
  * Density
  * Specific entropy
  * Specific internal energy
  * Specific isobaric heat capacity
  * Specific isochoric heat capacity
  * Speed of sound
  * Viscosity
  * Vapour fraction

All properties can be calculated with the inputs, p and T known, p and h known,
h and s known and some with pressure and density known.

XSteam is a full implementation of the IF-97 formulation including all regions
and all backward functions for good calculation speed. The code is speed
optimized with pressure and enthalpy as inputs for dynamic simulations.

Examples:

```matlab
>> XSteam('h_pt', 1, 20) # returns the enthalpy of water at 1 bar and 20 °C
>> 84.0118 kJ/kg

>> XSteam('rho_ph', 1, 3000) # returns the density of steam at 1 bar and 3000 kJ/kg
>> 0.4056 kg/m3

>> XSteam('w_pt', 1, 20) # returns the speed of sound at 1 bar and 20  °C
>> 1483.4 m/s

>> XSteam('tSat_p', 1) returns the saturation temperature at 1 bar.
>> 99.6059  °C
```

The X Steam Tables are the perfect tool both for replacing paper tables and for
advanced calculations. The XSteam tables are open source and free of charge.