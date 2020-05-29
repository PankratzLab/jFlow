jFlow
====================
Efficient Visualization, Processing, and Annotation of Flow Cytometry Data

----------------------
jFlow is a robust software package that provides a one-stop solution for visualizing, gating, batch-processing, and annotating flow cytometry files. 

Three main software packages are provided: jFlow, FCSProcessingPipeline, and FlowAnnotator:

Package        				|  description  
----------------------------|--------------  
jFlow						|  interactive visualization program for reviewing and editing gating   
FCS Processing Pipeline 	|  batch processing pipeline for creating images en-masse, or calculating counts and percentages  
Flow Annotator				|  program for reviewing the images created by FCSProcessingPipeline, allowing for fast and efficient review of relevant cell populations  


Running jFlow
=============


Main Class					|  description  
----------------------------|--------------  
jFlow						|  `org.genvisis.fcs.jFlow`   
FCS Processing Pipeline 	|  `org.genvisis.fcs.auto.FCSProcessingPipeline`    
Flow Annotator				|  `org.genvisis.flowannot.FlowAnnotator`   

 - *The requirements for each main class can be found by running each class with the -h flag.*
 - *Summaries and examples of required XML files can be found in the [FileDefs](https://github.com/PankratzLab/jFlow/blob/master/readme/FileDefs.md) document.*


Building with Maven
=====================

Genvisis uses [Maven](https://maven.apache.org/) to manage its dependencies, build process and more.

Genvisis is structured as a [multi-module](https://maven.apache.org/guides/introduction/introduction-to-the-pom.html#Project_Aggregation) Maven project with the following structure:

| pom-jflow  
| '-- Genvflow  

The key things to know about this structure are:

* The top-level `pom-jflow` pom.xml manages the build order of submodules and all dependency versions.
* The `PLab-public` project contains common utility classes used by jFlow and the main [Genvisis](https://github.com/PankratzLab/Genvisis) project.
* Each individual project is laid out according to the [standard Maven directory layout](https://maven.apache.org/guides/introduction/introduction-to-the-standard-directory-layout.html).

## Getting started

### Command line

For building on the command line, simply run `mvn` from the top-level `pom-jflow` directory. This will call the [install goal](https://maven.apache.org/guides/introduction/introduction-to-the-lifecycle.html) by default - building all modules and copying them to your [local Maven repository](https://maven.apache.org/guides/introduction/introduction-to-repositories.html).

### Eclipse

1. Import the top-level jFlow directory as an [Existing Maven Project](http://javapapers.com/java/import-maven-project-into-eclipse/). Eclipse will create projects automatically for `pom-jflow` and all sub-modules.

That's it! You can now develop jFlow code.

#### Caveats / Troubleshooting

When importing, Eclipse will recognize that jFlow is a multi-module build and create a project for each module. However, if additional modules are added later, Eclipse will not automatically import these.

So if you know a new module has been added, or are seeing odd missing dependency errors, try re-importing Maven projects from the `pom-jflow` directory.

#### Error: Project configuration is not up-to-date with pom.xml

Maven `pom.xml`s are translated to Eclipse projects via the M2Eclipse plugin. Unfortunately, by default it does not automatically update settings when the `pom.xml` changes, resulting in these annoying errors popping up from time to time.

While it is safe to simply update your project (e.g. with the [Quick Fix](http://help.eclipse.org/neon/index.jsp?topic=%2Forg.eclipse.jdt.doc.user%2Fconcepts%2Fconcept-quickfix-assist.htm) feature), there is also an Eclipse preference to [automatically update](http://www.eclipse.org/m2e/documentation/release-notes-16.html#new-experimental-auto-45-update-configuration-feature) - which can save some time and confusion.

#### Creating a run configuration

For convenience, you can also [create a dedicated run configuration](https://www.genuitec.com/products/myeclipse/learning-center/maven/launch-maven4myeclipse-maven-run-setup-tutorial/#2_Creating_a_CustomMavenLaunch_Configuration) for the top-level `pom-jflow` project.

Creating a run configuration will allow you to select the Maven build from the [Run toolbar button](https://developers.google.com/eclipse/docs/running_and_debugging_2_0), and also allows customization of properties (in contrast to the global properties declared in your `settings.xml`).

#### Style templates

This repository includes Eclipse style templates, located in `jFlow.git/config`. Before committing changes to jFlow, please import these templates to your [Clean Up](https://help.eclipse.org/neon/index.jsp?topic=%2Forg.eclipse.jdt.doc.user%2Freference%2Fpreferences%2Fjava%2Fcodestyle%2Fref-preferences-cleanup.htm) and [Formatter](https://help.eclipse.org/neon/index.jsp?topic=%2Forg.eclipse.jdt.doc.user%2Freference%2Fpreferences%2Fjava%2Fcodestyle%2Fref-preferences-formatter.htm) code style preferences.

To avoid the burden of remembering to manually format code, after setting the Clean Up and Formatter templates, you can also tell Eclipse to [automatically apply these formatters](https://stackoverflow.com/a/15655278/1027800) when saving files.

Synchronizing code styles across all developers makes reviewing code changes much easier, as the reader is not forced to try and separate formatting vs semantic changes.

In the event that the style templates themselves require updating, such a commit should also the new templates to the complete code base (ensuring a fresh starting point).

## Building the jFlow jar

Whether from Eclipse or the command line, to build the `genvflow.jar`, we need to run Maven with the parameters:

* project: **pom-jFlow** (or top directory of git repo)
* goal: **install**
* profile: **flow** 

When complete, your `genvflow.jar` will be built in the `Genv-flow/target/` directory, per the standard directory layout.

Notes:
- Without the `fastTests` profile enabled, the unit tests will be prohibitively slow when building.
- Without the `flow` profile enabled, the `genvflow.jar` will not be built (but this will speed up the Maven build time)

These parameters can be overridden [the standard Maven way](http://books.sonatype.com/mvnref-book/reference/running-sect-options.html) from the command line to customize the output .jar.

Additionally, when creating a run configuration in Eclipse there is a section to override parameters, so you can use multiple run configurations to manage your different views of the  codebase.

### Changing main class or jar name

The following maven properties can also be set when building the `.jar`:

* `outname`  - name of the .jar, default: "genvflow"
* `mainClass` - fully qualified main class to run when executing the .jar, default: org.genvisis.flowannot.FlowAnnotator

## Running tests

The default when a maven install is executed is to run every test. In order to have a faster build that skips lengthy tests, the `fastTests` profile can be activated by setting the `fastTests` property to `true`. (`-DfastTests=true` from the command line or added as a property in an eclipse run conifguration)

