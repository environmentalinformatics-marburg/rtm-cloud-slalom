## Important ##
The following instructions only refer to CLOUD and SLALOM versions published before Jun 23, 2011. From this date on, CLOUD and SLALOM are no longer implemented as two separate programs but combined into RTM\_CLOUD\_SLALOM. If you do not really want to use the deprecated versions, please download a version after Jun 23, 2011 and refer to CompilationAndRunning.

## Prerequisite ##

For running CLOUD or SLALOM you need
  * any Fortran 90 compiler
  * about 2,5 GB memory
  * about 2,5 GB of hard disk space.
<br></li></ul>


<h2>Downloading CLOUD, SLALOM and LUTs ##

CLOUD and SLALOM can be compiled and used independently of each other. Hence, if you are only interested in radiative transfer simulations, just go for CLOUD and if you are only interested in cloud property retrievals, go for SLALOM.

Both CLOUD and SLALOM require several look-up tables (LUTs). Therefore you have to download the programs and the look-up tables.

To get the program, please get the source code from
  * the [Download](http://code.google.com/p/rtm-cloud-slalom/downloads/list) page or
  * the [Source](http://code.google.com/p/rtm-cloud-slalom/source/checkout) page.

The Source page generally offers the latest versions of the code on the expense of potential bugginess. In contrast to the Download section, the source code of both programs is packed in one repository.

To get the look-up tables, please
  * got to the Download page and download all `rtm-cloud-slalom-lut.zip.001` to `rtm-cloud-slalom-lut.zip.004` files
  * unpack the `rtm-cloud-slalom-lut.zip.001` file to a folder of your choice (this will also unpack the other zip files).
<br></li></ul>


<h2>Compiling and Running CLOUD and SLALOM ##

If you get the program from the download section, unpack the zip file.

To compile the program, just navigate to the src folder and compile
  * `cloud.f90` and/or
  * `slalom.f90`

No compiler flags should be necessary for any compiler. For example, if you use gfortran just type:
```
gfortran -o slalom slalom.f90
```

<br>

<b>Running and Testing CLOUD</b>
To run CLOUD, you have to adapt the settings in the <code>cloud.cfg</code> file. This also includes the path to the LUTs which you should have downloaded already.<br>
<br>
To check your download and compilation, adapt the path to the LUTs in the  <code>cloud.cfg</code> file without any other changes and execute the program.<br>
<br>
Compare the content of the freshly generated <code>solar_lut_0856_10.dat</code> file with the <code>solar_lut_0856_10_test.dat</code> file supplied in the src folder.<br>
<br>
<br>

<b>Running and Testing SLALOM</b>
To run SLALOM, you do not need to adapt any configuration file but you have to execute the program with some command line parameters. To get a list of all necessary parameters, just run the program without any parameter.<br>
<br>
To check your download and compilation, execute SLALOM with the following parameters:<br>
<br>
<pre><code>./slalom 0 1 0.856 1.630 &lt;path-to-lut-folder&gt;/ lut_RInf_0856_aef_10.dat lut_RInf_1630_aef_10.dat 1 input_correct_20080520.dat output_correct_20080520.dat 1<br>
</code></pre>

Compare the content of the freshly generated <code>output_correct_20080520.dat</code> file with the <code>output_correct_20080520_test.dat</code> file supplied in the src folder.