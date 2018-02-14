#!/bin/sh
#This script is intended for launch on *nix machines
#-Xmx2000m indicates 2000 mb of memory, adjust number up or down as needed
#Script must be in the same directory as vis.jar
prefix=`dirname $(readlink $0 || echo $0)`
exec java -Xmx2000m \
	-cp "$prefix"/genvisis.jar cnv.Launch "$@"

