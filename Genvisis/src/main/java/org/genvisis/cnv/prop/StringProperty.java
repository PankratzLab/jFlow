package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.GROUP;

public class StringProperty extends Property<String> {
	public StringProperty(Project proj, String name, String desc, GROUP group, boolean editable,
												String defVal) {
		super(proj, name, desc, group, editable, defVal);
	}

	@Override
	public void parseValue(String valueStr) {
		setValue(valueStr);
	}
}
