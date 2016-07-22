Genvisis
====================
Efficient Visualization of Genomic data

----------------------
Genvisis is a robust software package that takes advantage of unique compressed data structures to efficiently process, assess quality of, analyze, and visualize the intensity data from GWAS arrays.

Specialized algorithms are used to call CNVs, and each algorithm uses its own approach. Genvisis allows a CNV to be inspected by visualizing boundaries overlaid with probe intensities


View the Change log [here](https://github.com/npankrat/Genvisis/blob/master/CHANGELOG.md)
 

[//]: # (1.	Parse Illumina/Infimax/csv/)
[//]: # (Load and parse micro array data (Create a new project))
[//]: # (2.	Extract Plots:)
[//]: # (Convert the data into a format that is)
[//]: # (3.	Slim Plots:)
[//]: # (An optional step to speed up the scatter plot for mid-size data)
[[//]: # (4.	Plots or Tools:)
[//]: # (Select the kind of analysis to run.)

[//]: # (While there are graphical software being so popular today, Genvisis provides some features that are not easily found elsewhere.)
[//]: # (Genvisis is a software product that Dr. Nathan Pankratz has developed for years. He started the project because he couldn’t find the right tool for his genetic research.)

Building with Maven
=====================

Genvisis uses [Maven](https://maven.apache.org/) to manage its dependencies, build process and more.

Genvisis is structured as a [multi-module](https://maven.apache.org/guides/introduction/introduction-to-the-pom.html#Project_Aggregation) Maven project with the following structure:

| pom-genvisis  
| '-- Logicle  
| '-- CFCS  
| '-- Genvisis  
| '-- Genv-experimental  
| '-- Genv-deprecated  
| '-- Assembly  

The key things to know about this structure are:

* The top-level `pom-genvisis` pom.xml manages the build order of submodules and all dependency versions.
* The `Assembly` project controls how the bundled `Genvisis.jar` is created. Any dependencies to this project will be included in the final .jar, and the pom.xml version number is used as the overall application version.
* Each individual project is laid out according to the [standard Maven directory layout](https://maven.apache.org/guides/introduction/introduction-to-the-standard-directory-layout.html).

## Getting started

### Command line

For building on the command line, simply run `mvn` from the top-level `pom-genvisis` level. This will call the [install goal](https://maven.apache.org/guides/introduction/introduction-to-the-lifecycle.html) by default - building all modules and copying them to your [local Maven repository](https://maven.apache.org/guides/introduction/introduction-to-repositories.html).

### Eclipse

1. Import the top-level Genvisis directory as an [Existing Maven Project](http://javapapers.com/java/import-maven-project-into-eclipse/). Eclipse will create projects automatically for `pom-genvisis` and all sub-modules.

That's it! You can now develop Genvisis code.

#### Building the Genvisis application

To build the `Genvisis.jar`, run either the `Assembly` or the `pom-genvisis` (if you need to rebuild other modules) project [as a Maven build](https://books.sonatype.com/m2eclipse-book/reference/running-sect-running-maven-builds.html), and select the **install** goal.

All Maven output will print in the `Console` tab. When complete, your `Genvisis.jar` will be built in the `Assembly/target/` directory, per the standard directory layout.

## Automatic uploads

You can use Maven to automatically upload `genvisis.jar` to a remote location after each build. To enable this feature, create a `settings.xml` in your `${user.home}/.m2` directory with the following structure:

```xml
<settings>
	<profiles>
		<profile>
			<id>genvisis-upload</id>
			<properties>
				<genv.upload.exe>[p]scp</genv.upload.exe>
				<genv.upload.path>user[:pass]@host:/path/to/output/</genv.upload.path>
			</properties>
		</profile>
	</profiles>

	<activeProfiles>
		<activeProfile>genvisis-upload</activeProfile>
	</activeProfiles>
</settings>
```

Once the `genv.upload.path` is set, every time the `Assembly` component builds the `Install` step, it will use this information to scp `genvisis.jar` to the specified remote path.
