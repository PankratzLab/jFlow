XML File Definitions  
======================

# JFlow

`jFlow` requires no XML files.

# Flow Annotator

### `panels=`

`FlowAnnotator` requires a `panels.xml` file as the first and only argument to the program  
  - This is the same file as required by the `FCSProcessingPipeline`.   
  - If this argument is not provided, the default included `panels.xml` file will be used, which defines two panels 

# FCS Processing Pipeline

`FCSProcessingPipeline` accepts a number of arguments, including:

### `pipe=`

Option     |  description  
-----------|-------------  
VIZ        |  Iterate through all linked FCS+WSP+additional-file-sets and generate visual plots for each cell population.  
BOOL	   |  Iterate through each linked file-set and generate inclusion matrices for (cell populations) * (data points)     
PCTS_CNTS  |  Iterate through each file-set and for each cell, compute the number of cells in each population as a raw integer and as the percentage of cells in the parent population  
		
### `panels=`

 - [default file](https://github.com/PankratzLab/jFlow/tree/master/Genv-flow/src/main/resources/docs/panels.xml)
 
Defines the structure of the wsp files for `FlowAnnotator` and `FCSProcessingPipeline`.  
Must follow the following specification:

```
<?xml version="1.0" encoding="UTF-8"?>
<panels>
	<!-- any number of panel tags may be present-->
	<panel>
		<!-- panel display name -->
		<name>Panel Example</name>
		<!-- any number of unique panel aliases, used to identify which FCS/WSP/image files belong to each panel -->
		<alias>panel_example</alias>
		<alias>panel example</alias>
		<alias>PanelExample</alias>
		<alias>pE</alias>
		<gateTree>
			<!-- the gating tree defines the structure of the image files created from the FCS/WSP files by the FCSProcessingPipeline. -->
			<!-- each self-contained node element defines a 'name' tag, and a 'parent' tag that matches a previous 'name' tag -->
			<!-- the first node is the root node, and doesn't define a parent -->
			<node name="boundary" />
			<node name="nonDebris" parent="boundary" />
			<node name="Lymphocytes (SSC-A v FSC-A)" parent="nonDebris" />
		</gateTree>
		<!-- each panel can define additional images to review, as composites of existing images -->
		<specials>
			<special>
				<!-- the key tag must match a 'name' tag from either the gateTree or addlFiles.xml file -->
				<key>Special Composite Name</key>
				<!-- value tags correspond to image names -->
				<value>effector memory helper Tcells (CCR7- CD45RA-)</value>
				<value>effector memory helper Tcells (CD95- CD28-)</value>
				<value>naive helper Tcells (CCR7/CD45RA)</value>
				<value>effector helper Tcells (CCR7/CD45RA)</value>
				<value>central memory helper Tcells (CCR7/CD45RA)</value>
				<value>effector memory helper Tcells (CCR7/CD45RA)</value>
			</special>
		</specials>
	</panel>
</panels>
```

### `dimensionOverrides=`

 - [default file](https://github.com/PankratzLab/jFlow/tree/master/Genv-flow/src/main/resources/docs/fluorochromeNameOverrides2.xml)  
 
Contains replacement names for fluorochrome channels.  Useful if data has been created from multiple machines with slightly altered fluorochrome names.

Example:

```
<?xml version="1.0" encoding="UTF-8"?>
<fluorochromes>
	<replace from="Comp-BV 605-A (CD95)" to="Comp-BV605-A (CD95)" />
    <replace from="Comp-BV 510-A (CD28)" to="Comp-BV510-A (CD28)" />
</fluorochromes>
```

This can be used to coalesce and combine names, e.g.:

```
    <replace from="Comp-BB 515-A (CD27)" to="Comp-BB515-A (CD27)" />
    <replace from="Comp-BB515-A (CD27)" to="Comp-FITC-A (CD27)" />
```


### `gateOverrideMatchup=`

 - [default file](https://github.com/PankratzLab/jFlow/tree/master/Genv-flow/src/main/resources/docs/overrideMatch.xml)

 Used to link cluster override files defined by the `gateOverrideDir` and `gateOverrideSfx` arguments to gate/population names in-memory.  
 
 Multiple `assign` tags act as a "union" operation, only allowing data points present in each of the given populations.  
   
 Example:
 
```
  
<?xml version="1.0" encoding="UTF-8"?>
<overrides>
	<override>
		<gate>Helper Tcells-CD4+</gate>
		<assign>HELPER_T</assign>
	</override>
	<override>
		<gate>naive helper Tcells (CD95- CD28+)</gate>
		<assign>HELPER_T</assign>
		<assign>naive</assign>
	</override>
</overrides> 
```

	
### `addlImgs=`

 - [default file](https://github.com/PankratzLab/jFlow/tree/master/Genv-flow/src/main/resources/docs/addlImgs.xml)  
 - Only used by the `VIZ` pipeline
 
Defines additional images to create, their names, parent populations, and x/y axes

Example:

```	
<?xml version="1.0" encoding="UTF-8"?>
<images>
	<image>
		<name>effector helper Tcells (CCR7/CD45RA)</name>
		<parent>effector helper Tcells (CCR7- CD45RA+)</parent>
		<x>Comp-BV 421-A (CCR7)</x>
		<y>Comp-BV 711-A (CD45RA)</y>
	</image>
</images>	
```

### `clusterOverrides=`  

 - [default file](https://github.com/PankratzLab/jFlow/tree/master/Genv-flow/src/main/resources/docs/clusterOverrides.xml)  
 - Only used by the `VIZ` pipeline
 
 Loads specified cluster assignments as data point colors.
 Example:
 
```
<?xml version="1.0" encoding="UTF-8"?>
<clusters>
	<cluster>
		<!-- 'parent' includes all cells in the parent gate; 'self' includes only cells in the specified gates -->
		<gate>parent</gate>
		<!-- gate/population names -->
		<key>effector memory helper Tcells (CCR7- CD45RA-)</key>
		<key>effector memory helper Tcells (CD95- CD28-)</key>
		rKey>naive helper Tcells (CCR7/CD45RA)</colorKey>
		<colorKey>effector helper Tcells (CCR7/CD45RA)</colorKey>
		<colorKey>central memory helper Tcells (CCR7/CD45RA)</colorKey>
		<colorKey>effector memory helper Tcells (CCR7/CD45RA</colorKey>
	</cluster>
</clusters>
```

