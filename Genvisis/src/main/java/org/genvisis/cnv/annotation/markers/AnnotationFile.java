package org.genvisis.cnv.annotation.markers;

import org.genvisis.cnv.filesys.AnnotationCollection;
import org.genvisis.cnv.filesys.Project;

/**
 * @author lane0212
 *
 *         This filetype may be useful when storing/accessing annotations that typically represent
 *         an entire project and are not updated frequently...i.e a more static version of
 *         {@link AnnotationCollection} <br>
 *         Going to be starting off with blast results...
 *
 */
public abstract class AnnotationFile {
	protected Project proj;
	protected String annotationFilename;
	protected Annotation[] annotations;
	protected AnalysisParams[] params;

	/**
	 * @param proj
	 * @param types used for extracting data types
	 * @param keys used for building/verifying the vcf header hash
	 * @param annotationFilename
	 */
	public AnnotationFile(Project proj, String annotationFilename) {
		super();
		this.proj = proj;
		this.annotationFilename = annotationFilename;
	}

	public Annotation[] getAnnotations() {
		return annotations;
	}

	public void setAnnotations(Annotation[] annotations) {
		this.annotations = annotations;
	}

	public AnalysisParams[] getParams() {
		return params;
	}

	public void setParams(AnalysisParams[] params) {
		this.params = params;
	}

}
