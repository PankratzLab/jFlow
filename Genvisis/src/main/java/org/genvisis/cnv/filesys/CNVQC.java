package org.genvisis.cnv.filesys;

import java.io.Serializable;

import org.genvisis.cnv.qc.CNVariantQC;
import org.genvisis.common.SerializedFiles;

public class CNVQC implements Serializable {

	public static final long serialVersionUID = 1L;
	private final CNVariantQC[] cnVariantQCs;

	public CNVQC(CNVariantQC[] cnVariantQCs) {
		this.cnVariantQCs = cnVariantQCs;
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static CNVQC load(String filename, boolean jar) {
		return (CNVQC) SerializedFiles.readSerial(filename, jar, true);
	}

	public CNVariantQC[] getCnVariantQCs() {
		return cnVariantQCs;
	}

}
