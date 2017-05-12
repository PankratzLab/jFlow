package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.COPY;
import org.genvisis.cnv.filesys.Project.GROUP;

public class BooleanProperty extends Property<Boolean> {
	public BooleanProperty(Project proj, String name, String desc, GROUP group, boolean editable,
												 COPY copyOnCorrection, Boolean defVal) {
		super(proj, name, desc, group, editable, copyOnCorrection, defVal);
	}

	@Override
	public void parseValue(String valueStr) {
		Boolean newValue = valueStr.equals("") ? getDefaultValue() : Boolean.valueOf(valueStr);
		setValue(newValue);
	}
}
