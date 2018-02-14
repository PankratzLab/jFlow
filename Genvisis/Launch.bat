#This script is intended for launch on Windows machines
#-Xmx2000m indicates 2000 mb of memory, adjust number up or down as needed
#Script must be in the same directory as vis.jar
for %%x in (%0) do set BatchPath=%%~dpsx
for %%x in (%BatchPath%) do set BatchPath=%%~dpsx
java  -Xmx22000m -cp %BatchPath%/genvisis.jar cnv.Launch  %*
PAUSE

