Scripts and pipelines that are specific to the Mouse Imaging Centre lab. These programs are not very likely to be applicable to other centres.

Contains:
- script to crop brains for the large field of view in-vivo scans
- script to convert CT images to MINC (adding scan information to the header)
- MEMRI / saddle coil reconstruction script
- distortion correction script (ex-vivo, saddle coil, Bruker cryo coil)

Prerequisites
-------------
Python3 

Installation
------------
The perl and C++ code:
<pre><code>
./autogen.sh
./configure
make
make install
</pre></code>

Python scripts:
<pre><code>
python3 python/setup.py install
</pre></code>


