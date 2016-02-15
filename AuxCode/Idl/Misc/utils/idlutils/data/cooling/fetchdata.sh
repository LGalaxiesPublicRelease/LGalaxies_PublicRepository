#!/bin/sh
# Download Dopita & Sutherland (1993) data
URL='http://www.mso.anu.edu.au/~ralph/data/cool'
curl --remote-name $URL/ABOUT.txt
curl --remote-name $URL/m+05.cie
curl --remote-name $URL/m-00.cie
curl --remote-name $URL/m-05.cie
curl --remote-name $URL/m-10.cie
curl --remote-name $URL/m-15.cie
curl --remote-name $URL/m-20.cie
curl --remote-name $URL/m-30.cie
curl --remote-name $URL/mzero.cie
echo Finished at `date`
