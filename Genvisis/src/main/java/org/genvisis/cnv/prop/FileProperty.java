package org.genvisis.cnv.prop;

import java.io.File;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.COPY;
import org.genvisis.cnv.filesys.Project.GROUP;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class FileProperty extends StringProperty {
	final boolean isDir;

	public FileProperty(Project proj, String name, String desc, GROUP group, boolean editable,
											COPY copyOnCorrection,
											String defVal, boolean dirOnly) {
		super(proj, name, desc, group, editable, copyOnCorrection,
					dirOnly ? ext.verifyDirFormat(defVal)
								 : ext.replaceAllWith(defVal, "\\",
																			"/")/*
																					 * == null || "".equals(defVal) ? null : new File(defVal)
																					 */);
		isDir = dirOnly;
	}

	public boolean isDirectory() {
		return isDir;
	}

	@Override
	public void setValue(String value) {
		super.setValue(isDir ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/"));
	}

	@Override
	public String getValueString() {
		return getValue();
	}

	@Override
	public String getValue() {
		return getValue(false, false);
	}

	public String getValue(boolean mkdirs, boolean verbose) {
		return getValue(null, mkdirs, verbose);
	}

	// TODO NP asks: When is this subdir option ever used?
	// TODO use ext.isrelativepath, document that this is a relative path, and ensure it's being used
	// as such
	public String getValue(String subdir, boolean mkdirs, boolean verbose) {
		String valu = super.getValue();
		if (valu.contains("~")) {
			valu = ext.replaceTilde(valu);
		}

		String tempValue = valu;
		// if (!"".equals(valu) && !valu.startsWith(".") && !valu.startsWith("/") && valu.indexOf(":")
		// == -1) {
		if (!"".equals(valu) && Files.isRelativePath(valu)) {
			if (isDir) {
				if (getName().equals(PropertyKeys.KEY_PROJECT_DIRECTORY)) {
					tempValue = valu;
					// jlanej removed the below since most MitoPipeline examples use relative paths, and
					// doesn't
					// play nice with $USER/.genvisis , can focus on
					// example properties more specifically later
					// LaunchProperties.directoryOfLaunchProperties() +
					// throw new IllegalArgumentException();


				} else {
					tempValue = getProject().PROJECT_DIRECTORY.getValue() + valu;
				}
			} else {
				tempValue = getProject().PROJECT_DIRECTORY.getValue()
										+ (subdir == null ? "" : getProject().getProperty(subdir).getValueString())
										+ valu;
			}
		}
		if (!Files.exists(tempValue, getProject().JAR_STATUS.getValue())) {
			if (mkdirs/* && getProject().JAR_STATUS.getValue() */) {
				if (isDir) {
					(new File(tempValue)).mkdirs();
				} else {
					(new File(ext.parseDirectoryOfFile(tempValue))).mkdirs();
				}
			} else if (verbose) {
				if (isDir) {
					getProject().getLog().reportError("Error - directory '" + tempValue + "' does not exist");
				} else {
					if (!Files.exists(ext.parseDirectoryOfFile(tempValue),
														getProject().JAR_STATUS.getValue())) {
						getProject().getLog()
												.reportError("Error - the directory ('"
																		 + ext.parseDirectoryOfFile(tempValue)
																		 + "') of the file you're trying to access/create ('"
																		 + ext.removeDirectoryInfo(tempValue) + "') does not exist");
					} else {
						getProject().getLog().reportError("Error - file '" + tempValue + "' does not exist");
					}
				}
			}
		}

		return tempValue;
	}
}
