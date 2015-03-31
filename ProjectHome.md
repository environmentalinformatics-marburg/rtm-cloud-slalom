## Radiative Transfer Model CLOUD ##
CLOUD is a radiative transfer model for the simulation of reflection, transmission, and absorption characteristics of terrestrial clouds. It is fast, accurate, and capable of calculating multiple radiative characteristics of cloudy media including the spherical and plane albedo, reflection and transmission functions, absorptance as well as global and diffuse transmittance. The approximations are based on the asymptotic solutions of the radiative transfer equations valid at cloud optical thicknesses larger than 5.

<br>

<h2>Cloud Property Retrieval SLALOM</h2>
SLALOM is a daytime retrieval of cloud optical and microphysical parameters from optical satellite data. It is capable of retrieving the cloud optical thickness, the effective cloud droplet radius, the liquid and ice water paths, the particle absorption length as well as some other properties of water and ice clouds. The technique is based on simple yet highly accurate approximations of the asymptotic solutions of the radiative transfer theory which have already been implemented in the forward radiative transfer model CLOUD. These approximations enable a solution of the equations during run time leading to a very fast computation speed. Since these asymptotic solutions are generally applicable to weakly absorbing media only, pre-calculated look-up tables for the reflection function of a semi-infinite cloud (and also the escape function) are used to overcome this restriction within this new retrieval.<br>
<br>
<br>

<h2>Implementation</h2>
CLOUD and SLALOM are implemented in Fortran. Although CLOUD is an integral part of SLALOM, both programs have initially been developed and published as separate sources (see versions <= 1.0.1 in the download section). Starting from Jun 23, 2011, the two software packages have been combined into the RTM_CLOUD_SLALOM program. Hence, if one wants to use CLOUD or SLALOM can be defined during run time.<br>
<br>
For compiling and installation instructions please refer to the CompilationAndRunning page of our Wiki. For older revisions of the program with separate CLOUD and SLALOM routines, see CompilationAndRunningSeparate.<br>
<br>
<br>

<h2>References</h2>
The CLOUD model is described in Kokhanovsky, A. A. and Nauss, T.: Reflection and transmission of solar light by clouds: asymptotic theory, Atmospheric Chemistry and Physics 6, 5537-5545, doi:10.5194/acp-6-5537-2006, 2006 (see <a href='http://www.atmos-chem-phys.org/6/5537/2006/acp-6-5537-2006.html'>http://www.atmos-chem-phys.org/6/5537/2006/acp-6-5537-2006.html</a>).<br>
<br>
The SLALOM retrieval is described in Nauss, T. and Kokhanovsky, A. A.: Retrieval of warm cloud optical properties using simple approximations, Remote Sensing of Environment 115/6, 1317-1325, doi: 10.1016/j.rse.2011.01.010, 2011 (see <a href='http://www.sciencedirect.com/science/article/pii/S0034425711000307'>http://www.sciencedirect.com/science/article/pii/S0034425711000307</a>).<br>
<br>
<br>

<h2>Other stuff</h2>
Please also visit our station and remote sensing data processing framework at <a href='http://code.google.com/p/julendat/'>julendat</a>.