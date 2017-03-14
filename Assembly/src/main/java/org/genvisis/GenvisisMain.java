package org.genvisis;

import org.genvisis.cnv.Launch;
import org.genvisis.common.LauncherManifest;

/**
 * This is a dummy entry point that currently simply delegates to {@link Launch}. It was primarily
 * added because https://javafx-maven-plugin.github.io/ was unhappy about building a jar with no
 * source files (which, to be fair, the .jar assembly had been complaining about as well). This
 * could serve as a container for JavaFX-specific launch considerations.
 */
public class GenvisisMain {

	public static void main(String... args) {
		LauncherManifest.setLaunchClass(GenvisisMain.class);
		Launch.main(args);
	}
}
