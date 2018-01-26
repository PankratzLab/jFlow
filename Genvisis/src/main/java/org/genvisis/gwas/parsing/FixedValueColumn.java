package org.genvisis.gwas.parsing;

import java.util.Map;

public class FixedValueColumn extends AbstractFileColumn<String> {

	private String value;

	public FixedValueColumn(String name, String value) {
		super(name);
		this.value = value;
	}

	@Override
	public String getValue(String[] line) throws ParseFailureException {
		return value;
	}

	@Override
	public void initialize(Map<String, Integer> headerMap) {
		// no-op
	}


	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((getName() == null) ? 0 : getName().hashCode());
		result = prime * result + value.hashCode();
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		FixedValueColumn other = (FixedValueColumn) obj;
		if (!value.equals(other.value))
			return false;
		if (getName() == null) {
			if (other.getName() != null)
				return false;
		} else if (!getName().equals(other.getName()))
			return false;
		return true;
	}

}
