## Important ##
The following instructions only refer to software versions after Jun 23, 2011, where CLOUD and SLALOM are combined into the new RTM\_CLOUD\_SLALOM program.
If you are interested in software versions before Jun 23, 2011, please refer to CompilationAndRunningSeparate.

## Prerequisite ##

For running RTM\_CLOUD\_SLALOM you need
  * any Fortran 90 compiler
  * about 2,5 GB memory
  * about 2,5 GB of hard disk space.
<br></li></ul>


<h2>Downloading RTM\_CLOUD\_SLALOM and LUTs ##

RTM\_CLOUD\_SLALOM combines the previously separate CLOUD and SLALOM program into one source code.

The program requires several look-up tables (LUTs). Therefore you have to download the programs and the look-up tables.

To get the program, please get the source code from
  * the [Download](http://code.google.com/p/rtm-cloud-slalom/downloads/list) page or
  * the [Source](http://code.google.com/p/rtm-cloud-slalom/source/checkout) page.

The Source page generally offers the latest versions of the code on the expense of potential bugginess.

To get the look-up tables, please
  * got to the Download page and download all `rtm-cloud-slalom-lut.zip.001` to `rtm-cloud-slalom-lut.zip.004` files
  * unpack the `rtm-cloud-slalom-lut.zip.001` file to a folder of your choice (this will also unpack the other zip files).
<br></li></ul>


<h2>Compiling and Running RTM\_CLOUD\_SLALOM ##

If you get the program from the download section, unpack the zip file.

To compile the program, just navigate to the src folder and compile
`rtm_cloud_slalom.f90`

No compiler flags should be necessary for any compiler. For example, if you use gfortran just type:
```
gfortran -o rtm_cloud_slalom rtm_cloud_slalom.f90
```

<br>

<b>Running and Testing the forward model CLOUD</b>
To run CLOUD, you have to adapt the settings in the <code>cloud.cfg</code> file. This also includes the path to the LUTs which you should have downloaded already.<br>
<br>
To check your download and compilation, adapt the path to the LUTs in the  <code>cloud.cfg</code> file without any other changes and execute RTM_CLOUD_SLALOM in CLOUD mode with the following parameter:<br>
<pre><code>./rtm_cloud_slalom cloud<br>
</code></pre>


Compare the content of the freshly generated <code>solar_lut_0856_10.dat</code> file with the <code>solar_lut_0856_10_test.dat</code> file supplied in the src folder.<br>
<br>
<br>

<b>Running and Testing the cloud property retrieval SLALOM</b>
To run SLALOM, you do not need to adapt any configuration file but you have to execute the program with some command line parameters. To get a list of all necessary parameters, just run the program without any parameter.<br>
<br>
To check your download and compilation, execute RTM_CLOUD_SLALOM in SLALOM mode with the following parameters:<br>
<br>
<pre><code>./rtm_cloud_slalom slalom 0 1 0.856 1.630 &lt;path-to-lut-folder&gt;/ lut_RInf_0856_aef_10.dat lut_RInf_1630_aef_10.dat 1 input_correct_20080520.dat output_correct_20080520.dat 1<br>
</code></pre>

Compare the content of the freshly generated <code>output_correct_20080520.dat</code> file with the <code>output_correct_20080520_test.dat</code> file supplied in the src folder.