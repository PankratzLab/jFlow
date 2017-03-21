package org.genvisis.cnv.annotation.markers;

import org.genvisis.cnv.filesys.AnnotationCollection;
import org.genvisis.common.Logger;

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
	protected String annotationFilename;
	protected Annotation[] annotations;
	protected AnalysisParams[] params;
	protected Logger log;

	/**
	 * @param annotationFilename
	 * @param log 
	 */
	public AnnotationFile(String annotationFilename, Logger log) {
		super();
		this.annotationFilename = annotationFilename;
		this.log = log;
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
