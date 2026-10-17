# Flux Spreader

Shortwave radiative fluxes are spread across fine timesteps to calculate mean values over longer periods (days or months). This approach consists in a virtual rescaling of the pre-computed fluxes to an evolving distrubition of the solar irradiance, without the need to recompute fluxes at each timestep. In fact, even in models such as the IFS, the radiation scheme is called only every few timesteps (typically every hour, while the model timestep might be on the order of ten minutes) in order to keep the computational cost of radiation calculations manageable.

When performing offline calculations, the coarse temporal resolution is often imposed by the availability of model data (e.g., hourly or three-hourly outputs) or by the need to reduce the total amount of data handled and calculations performed, and spreading fluxes across finer timesteps allows nonetheless an effective calculation of daily- or monthly-mean fluxes.

## Solar angles
The solar zenith angle $\theta$ is the angle between the local vertical and the direction of the sun and is expressed as a function of the time of day, the day of the year, and the geographic location on the Earth. The cosine of the solar zenith angle, $\mu_0 = \cos \theta$, is used in radiative transfer calculations to account for the angle of incidence of solar radiation. For the calculation of reflectances seen by a sensor located in orbit, also the solar azimuth angle is needed. These astronomical quantities are calculated natively by ecRadOff within the eot (Equation of Time) module, which consists of a minimal class to support efficient calculation of solar angles (see API for {doc}`../api/eot`).

## Concept
The solar irradiance can be expressed as a function of time, latitude ($\phi$), and longitude ($\lambda$):

$I(t, \phi, \lambda) = I_0 \mu_0(t, \phi, \lambda)$

where $I_0$ is the solar constant and $\theta$ is the solar zenith angle at any given time and location. Normally, $\mu_0 = \cos \theta \in [0, 1]$ is referred to as the cosine of the solar zenith angle. Replacing $\mu_0$  with a smooth, strictly non-zero function $\mu_0^*$ allows the calculation of non-zero fluxes even when the sun is below the horizon, that can eventually be re-scaled to the updated $\mu_0(t', \phi, \lambda)$ at a subsequent timestep $t'$. Let the flux at time $t$ $F(t, \phi, \lambda)$, then:

$F(t', \phi, \lambda) \approx F(t, \phi, \lambda) \frac{\mu_0(t', \phi, \lambda)}{\mu_0^*(t, \phi, \lambda)}$

$\mu_0^*$ can be calculated by calling :func:`eot.cos_sza_rad` with the parameter ``night_to_zero=False``.

## Usage
Support for flux spreading is implemented...