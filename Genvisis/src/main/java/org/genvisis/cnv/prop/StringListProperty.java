package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.GROUP;
import org.genvisis.common.Array;
import org.genvisis.common.ext;

public class StringListProperty extends Property<String[]> {
	static String delim = ";";
	boolean isFile = false;
	boolean isDir = false;

	public boolean isFile() {
		return isFile;
	}

	public boolean isDirectory() {
		return isDir;
	}

	public StringListProperty(Project proj, String name, String desc, GROUP group, boolean editable,
														String[] defVal, boolean file, boolean dir) {
		super(proj, name, desc, group, editable, defVal);
		isFile = file;
		isDir = dir;
	}

	public StringListProperty(Project proj, String name, String desc, GROUP group, boolean editable,
														String defVal, boolean file, boolean dir) {
		super(proj, name, desc, group, editable,
					defVal.equals("") ? new String[0] : defVal.split(delim));
		isFile = file;
		isDir = dir;
	}

	@Override
	public void parseValue(String valueStr) {
		String[] pts = "".equals(valueStr) ? new String[0] : valueStr.split(delim);
		setValue(pts);
	}

	@Override
	public String getValueString() {
		return Array.toStr(getValue(), delim);
	}

	@Override
	public String getDefaultValueString() {
		return Array.toStr(getDefaultValue(), delim);
	}

	public String getDelimiter() {
		return delim;
	}

	@Override
	public String[] getValue() {
		String[] values = super.getValue();
		if (isFile || isDir) {
			for (int i = 0; i < values.length; i++) {
				if (!"".equals(values[i])	&& !values[i].startsWith(".") && !values[i].startsWith("/")
						&& values[i].indexOf(":") == -1) {
					values[i] = getProject().PROJECT_DIRECTORY.getValue() + values[i];
				}
			}
		}
		return values;
	}

	@Override
	public void setValue(String[] value) {
		String[] values = value;
		if (isFile || isDir) {
			for (int i = 0; i < values.length; i++) {
				if (values[i].startsWith(getProject().PROJECT_DIRECTORY.getValue())) {
					values[i] = values[i].substring(getProject().PROJECT_DIRECTORY.getValue().length());
				}
			}
		}
		super.setValue(values);
	}

	public void removeValue(String value) {
		String[] newValues = new String[getValue().length - 1];
		String[] values = getValue();
		int index = 0;
		for (String value2 : values) {
			boolean skip = false;
			if ((isFile || isDir) && (ext.verifyDirFormat(value2).equals(ext.verifyDirFormat(value)))) {
				skip = true;
			} else if (value2.equals(value)) {
				skip = true;
			}
			if (skip) {
				continue;
			}
			newValues[index++] = value2;
		}
		setValue(newValues);
	}

	public void addValue(String value) {
		addValue(value, 0);
	}

	public void addValue(String valu, int index) {
		String[] curValue = getValue();

		for (String s : curValue) {
			// Check if the current value is already present
			if (s.equals(valu)) {
				return;
			}
		}

		String value = valu;
		if (isDir) {
			value = ext.verifyDirFormat(value);
		} else if (isFile) {
			value = ext.verifyDirFormat(value);
			value = value.substring(0, value.length() - 1);
		}
		setValue(Array.addStrToArray(value, curValue, index));
	}
}
