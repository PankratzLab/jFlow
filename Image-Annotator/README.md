## Building Image-Annotator distribution

1. Run `mvn clean install -PfastTests` from the root Genvisis directory to install the current Genvisis libraries in your local maven repository.
1. Run `mvn clean install -Prelease,jfx-distribution` from the Image-Annotator directory to build the Image-Annotator distribution.
1. Delete the current contents of `genvisis.org/dist-image-annotator/`
1. Upload the contents of `Image-Annotator/target/app/` to `genvisis.org/dist-image-annotator/`

The next time the Image-Annotator app is launched it will automatically download the changes.

## Building Image-Annotator installer

Building the installer is the same as building the distribution, except you add one more profile when building Image-Annotator:

`mvn clean install -Prelease,jfx-distribution,jfx-installer`

This will put an installer.exe in your `Image-Annotator/target/fxinstaller/bundle/` directory.

Currently this file is uploaded as `genvisis.org/image-annotator/image-annotator-installer.exe`

NOTE: the installer does not need to be rebuilt each time. It is simply a platform-specific launcher that knows where to look for updates. You would only need to rebuild it if you changed the app directory serving the updates.

NOTE: the installer does need to be built once for each version of Java to be distributed. It uses the version of the JDK used to build it.

## From Eclipse

The steps are the same from Eclipse as the command line. Note:

* Image-Annotator has to be imported into Eclipse manually, since it is not a child project of pom-genvisis.
* See these [instructions for running Maven builds from Eclipse](https://books.sonatype.com/m2eclipse-book/reference/running-sect-running-maven-builds.html)

You only need to set up these builds once in Eclipse, and then can select them from the list of existing run profiles.
