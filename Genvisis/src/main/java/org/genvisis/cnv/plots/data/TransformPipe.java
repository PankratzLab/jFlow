package org.genvisis.cnv.plots.data;


public abstract class TransformPipe extends AbstractPipe {

	public abstract String transformValue(String value);

	@Override
	public String pipeValue(String value) throws RejectedValueException {
		return transformValue(value); // Transforms don't reject values
	}
}
