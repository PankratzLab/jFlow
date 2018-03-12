Genvisis
====================
Efficient Visualization of Genomic data

----------------------
Genvisis is a robust software package that takes advantage of unique compressed data structures to efficiently process, assess quality of, analyze, and visualize the intensity data from GWAS arrays.

Specialized algorithms are used to call CNVs, and each algorithm uses its own approach. Genvisis allows a CNV to be inspected by visualizing boundaries overlaid with probe intensities.

Get started with the [Genvisis manual](https://github.com/PankratzLab/Genvisis/wiki/About-Genvisis).

View the Change log [here](https://github.com/npankrat/Genvisis/blob/master/CHANGELOG.md).
 

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

For building on the command line, simply run `mvn` from the top-level `pom-genvisis` directory. This will call the [install goal](https://maven.apache.org/guides/introduction/introduction-to-the-lifecycle.html) by default - building all modules and copying them to your [local Maven repository](https://maven.apache.org/guides/introduction/introduction-to-repositories.html).

### Eclipse

1. Import the top-level Genvisis directory as an [Existing Maven Project](http://javapapers.com/java/import-maven-project-into-eclipse/). Eclipse will create projects automatically for `pom-genvisis` and all sub-modules.

That's it! You can now develop Genvisis code.

#### Caveats / Troubleshooting

When importing, Eclipse will recognize that Genvisis is a multi-module build and create a project for each module. However, if additional modules are added later, Eclipse will not automatically import these.

So if you know a new module has been added, or are seeing odd missing dependency errors, try re-importing Maven projects from the `pom-genvisis` directory.

#### Building the Genvisis application

Just like on the command line, to build the `genvisis.jar`, run the `pom-genvisis` project [as a Maven build](https://books.sonatype.com/m2eclipse-book/reference/running-sect-running-maven-builds.html), and select the **install** goal.

All Maven output will print in the `Console` tab. When complete, your `genvisis.jar` will be built in the `Assembly/target/` directory, per the standard directory layout.

#### Error: Project configuration is not up-to-date with pom.xml

Maven `pom.xml`s are translated to Eclipse projects via the M2Eclipse plugin. Unfortunately, by default it does not automatically update settings when the `pom.xml` changes, resulting in these annoying errors popping up from time to time.

While it is safe to simply update your project (e.g. with the [Quick Fix](http://help.eclipse.org/neon/index.jsp?topic=%2Forg.eclipse.jdt.doc.user%2Fconcepts%2Fconcept-quickfix-assist.htm) feature), there is also an Eclipse preference to [automatically update](http://www.eclipse.org/m2e/documentation/release-notes-16.html#new-experimental-auto-45-update-configuration-feature) - which can save some time and confusion.

#### Creating a run configuration

For convenience, you can also [create a dedicated run configuration](https://www.genuitec.com/products/myeclipse/learning-center/maven/launch-maven4myeclipse-maven-run-setup-tutorial/#2_Creating_a_CustomMavenLaunch_Configuration) for the top-level `pom-genvisis` project.

Creating a run configuration will allow you to select the Maven build from the [Run toolbar button](https://developers.google.com/eclipse/docs/running_and_debugging_2_0), and also allows customization of properties (in contrast to the global properties declared in your `settings.xml`).

#### Changing main class or jar name

The Assembly project uses two parameters to determine how to build the final `.jar`:

* `outname`  - name of the .jar, default: "genvisis"
* `mainClass` - fully qualified main class to run when executing the .jar, default: org.genvisis.cnv.Launch

These parameters can be overridden [the standard Maven way](http://books.sonatype.com/mvnref-book/reference/running-sect-options.html) to customize the output .jar.

Additionally, when creating a run configuration in Eclipse there is a section to override parameters, so you can use multiple run configurations to manage your different views of the Genvisis codebase.

A final option is to create a `pom.xml` that uses the `Assembly` pom as its parent and overrides these parameters. The advantage to this method is that it provides a tangible `pom.xml` that implicitly documents the parameters, and can be shared with the community

#### Building a unified jar

Performing application assembly (merging all classes into a single unified jar) is costly, therefore this behavior must be explicitly enabled by turning profiles on. For example, the `genv` profile turns on 

#### Running tests

The default when a maven install is executed is to run every test. In order to have a faster build that skips lengthy tests, the `fastTests` profile can be activated by setting the `fastTests` property to `true`. (`-DfastTests=true` from the command line or added as a property in an eclipse run conifguration)

#### Style templates

This repository includes Eclipse style templates, located in `Genvisis.git/config`. Before committing changes to Genvisis, please import these templates to your [Clean Up](https://help.eclipse.org/neon/index.jsp?topic=%2Forg.eclipse.jdt.doc.user%2Freference%2Fpreferences%2Fjava%2Fcodestyle%2Fref-preferences-cleanup.htm) and [Formatter](https://help.eclipse.org/neon/index.jsp?topic=%2Forg.eclipse.jdt.doc.user%2Freference%2Fpreferences%2Fjava%2Fcodestyle%2Fref-preferences-formatter.htm) code style preferences.

To avoid the burden of remembering to manually format code, after setting the Clean Up and Formatter templates, you can also tell Eclipse to [automatically apply these formatters](https://stackoverflow.com/a/15655278/1027800) when saving files.

Synchronizing code styles across all developers makes reviewing code changes much easier, as the reader is not forced to try and separate formatting vs semantic changes.

In the event that the style templates themselves require updating, such a commit should also the new templates to the complete code base (ensuring a fresh starting point).

## Automatic upload and copy

The following properties are available when the `genv` profile is activated:

* `genv.copy.path` - the output `genvisis.jar` will be copied to this location
* `genv.upload.exe` - this application (e.g. `scp`, `pscp`) will be used to upload the output `genvisis.jar`
* `genv.upload.path` - the remote location to send the `genvisis.jar`, e.g. `user@host:/path/`

### Note on passwords

You should not store your password in plaintext in this settings file. The best option is to [set up a public/private key pair](https://help.ubuntu.com/community/SSH/OpenSSH/Keys) with the server you are uploading to. Another option is encryption with Maven directly via [settings-security.xml](https://maven.apache.org/guides/mini/guide-encryption.html)

## Creating a native application

Genvisis uses the [JavaFX-maven-plugin](https://github.com/javafx-maven-plugin/javafx-maven-plugin) to create native bundles. This configuration is stored in the `Assembly` project's `pom.xml`, in the `jfx` profile. It can be activated by running:

```bash
mvn clean package -P jfx
```

from either the `/Assembly/` directory or the project base directory (containing Genvisis' top-level `pom.xml`).

This will create several outputs of interest:

Genvisis.git</br>
|&nbsp;&nbsp;'-- _Assembly_</br>
|&nbsp;&nbsp;&nbsp;&nbsp;'-- _target_</br>
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'-- **genvisis.tar.gz** (compressed native app w/ runtime)</br>
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'-- _jfx_</br>
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'-- **app** (folder with jar + dependencies)</br>
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'-- _native_</br>
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'-- **genvisis** (folder with app + platform-specific runtime)</br>
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'-- **genvisis-xxx.[rpm|deb|exe|msi|app]** (platform-specific installer)</br>


Some notes:

* Native bundles will only be built for the local OS
* Each OS requires its own build tools to be installed (e.g. debian on Ubuntu, rpm on Linux, [InnoSetup](http://www.jrsoftware.org/isinfo.php) for `.exe` or [WiX](http://wixtoolset.org/) for `.msi` on Windows)
* We recommend reading [the JavaFX tutorials](http://docs.oracle.com/javafx/2/deployment/self-contained-packaging.htm) for more information about deployment options.

