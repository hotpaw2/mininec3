# mininec3 for Chipmunk Basic
## for radio antenna modeling (computational electrodynamics)

This repository contains a lightly modifed version
of the MININEC3.BAS program.
The MININEC3.BAS BASIC program was
originally developed by J.C.Logan and J.W.Rockway 
at Navel Ocean Systems Center in San Diego, CA.
MININEC3 is a method-of-moments computer program
for the analysis of thin-wire antennas.

This version of mininec3.bas has been minimally modified,
just enough to be compatible with Chipmunk Basic,
command-line version 3.6.8 (1.368.x) or later.
Chipmunk Basic is a Basic interpreter
for macOS and Linux (Raspberry Pi, et.al.).  
(mininec3.bas might also run under HotPaw Basic on the Apple iPad,
but extensive text console input and file management
on the iPad is a bit problematic.)
A few variable names, keywords, operators,
and file I/O statements had to be modified, due to minor
differences between QBASIC and Chipmunk Basic.
No GUI or graphics have been added with this modification.
The mininec3 application is still strictly command-line driven.
But since Chipmunk Basic on recent MacBooks
interprets Basic roughly 100X faster than
compiled QBASIC on a 8086 IBM PC, 
the maximun number of antenna wire segments
has been increased from 50 to 250.
This still allows the method-of-moments linear equation solver
to run in reasonable time, even in a purely interpreded Basic.

Information on the original MININEC3 program can be found at these web sites:

  https://ieeexplore.ieee.org/document/14088

  https://www.arrl.org/files/file/Technology/tis/info/pdf/9102018.pdf

  https://en.wikipedia.org/wiki/Numerical_Electromagnetics_Code#MININEC

  http://wb0dgf.com/mininec.htm

  http://www.vectorbd.com/bfd/antenna/mininec3.inf

  https://apps.dtic.mil/dtic/tr/fulltext/u2/a181682.pdf

MININEC3.BAS is U.S. Fedrally funded code reportedly now in the Public Domain.

Also included is a 3-line text file "mininec.inp",
modified from fixed column position numeric fields to CSV format.
If the mininec2 program finds a file with this name,
it uses the data in this file
for wire segment input, instead of inputting from the console.
Try this file with:
+1 for Free Space,
7.150 MHz,
1 source of: 5,100,0 ,
0 loads

Modified by:  rhn@nicholson.com  N6YWU  2019-10-20

Chipmunk Basic Info: http://www.nicholson.com/rhn/basic/

--
