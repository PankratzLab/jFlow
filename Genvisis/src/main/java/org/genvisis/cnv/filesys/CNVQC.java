package org.genvisis.cnv.filesys;

import java.io.Serializable;

import org.genvisis.cnv.qc.CNVariantQC;
import org.genvisis.common.Files;

public class CNVQC implements Serializable {

	public static final long serialVersionUID = 1L;
	private CNVariantQC[] cnVariantQCs;

	public CNVQC(CNVariantQC[] cnVariantQCs) {
		this.cnVariantQCs = cnVariantQCs;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static CNVQC load(String filename, boolean jar) {
		return (CNVQC) Files.readSerial(filename, jar, true);
	}

	public CNVariantQC[] getCnVariantQCs() {
		return cnVariantQCs;
	}

}