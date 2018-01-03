package org.genvisis.jfcs;

import java.io.UnsupportedEncodingException;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.Map;

public class FCSKeywords {
	byte[] textBytes;
	Map<String, String> raw;
	private String[] parameterNames;
	private String[] parameterNamesLong;

	protected FCSKeywords(byte[] textBytes) throws UnsupportedEncodingException {
		this.textBytes = textBytes;
		parseRaw(textBytes);
		parseParamNames();
	}

	// specialty methods:
	public int getParameterCount() {
		return getKeywordInt(FCS_KEYWORD.PAR.keyword);
	}

	public String getParameterLongName(int index) {
		return parameterNamesLong[index];
	}

	public String getParameterShortName(int index) {
		return parameterNames[index];
	}

	public double[] getParamScaling() {
		double[] scal = new double[getParameterCount()];
		for (int i = 0, count = getParameterCount(); i < count; i++) {
			scal[i] = 1;
			try {
				scal[i] = Double.parseDouble(getKeyword(FCS_KEYWORD.PnG.getKey(i + 1)));
			} catch (Exception e) {
			}
		}
		return scal;
	}

	public int getEventCount() {
		return getKeywordInt(FCS_KEYWORD.TOT.keyword)
					 - (hasKeyword(FCS_KEYWORD.ABRT.keyword) ? getKeywordInt(FCS_KEYWORD.ABRT.keyword) : 0)
					 - (hasKeyword(FCS_KEYWORD.LOST.keyword) ? getKeywordInt(FCS_KEYWORD.LOST.keyword) : 0);
	}

	public ByteOrder getEndianness() {
		String ord = raw.get(FCS_KEYWORD.BYTEORD.keyword);
		if ("4,3,2,1".equals(ord)) {
			return ByteOrder.BIG_ENDIAN;
		} else if ("1,2,3,4".equals(ord)) {
			return ByteOrder.LITTLE_ENDIAN;
		}
		return null;
	}

	public String getKeyword(FCS_KEYWORD keyword) {
		return getKeyword(keyword.keyword);
	}

	public String getKeyword(String keyword) {
		// TODO error checking / string checking (has $?)
		return raw.get(keyword);
	}

	public boolean hasKeyword(FCS_KEYWORD keyword) {
		return hasKeyword(keyword.keyword);
	}

	public boolean hasKeyword(String keyword) {
		// TODO error checking / string checking (has $?)
		return raw.containsKey(keyword);
	}

	public int getKeywordInt(String keyword) {
		// TODO lots of potential for error checking
		return Integer.parseInt(raw.get(keyword));
	}

	public long getKeywordLong(String keyword) {
		// TODO lots of potential for error checking
		return Long.parseLong(raw.get(keyword));
	}

	private void parseRaw(byte[] textBytes) throws UnsupportedEncodingException {
		byte delimByte = textBytes[0];

		raw = new HashMap<>();

		String v = null;
		byte[] sub;
		boolean isKey = true;
		int s = 1, e = 1;

		while (s < textBytes.length) {
			while (true) {
				e++;
				if (textBytes[e] == delimByte) {
					if (e == textBytes.length - 1 || textBytes[e + 1] != delimByte) {
						// delims can exist in values if doubled
						break;
					}
				}
			}
			sub = new byte[e - s];
			System.arraycopy(textBytes, s, sub, 0, e - s);
			if (!isKey) {
				// values can be UTF-8
				raw.put(v, new String(sub, isKey ? "US-ASCII" : "UTF-8").trim());
			} else {
				// keys are ascii
				v = new String(sub, "US-ASCII").trim();
				// trim() technically /shouldn't/ be necessary according to spec
			}
			isKey = !isKey;
			s = e + 1;
			e = e + 1;
		}
	}

	private void parseParamNames() {

		Map<String, String> paramNames = new HashMap<>();
		for (String k : raw.keySet()) {
			FCS_KEYWORD v = FCS_KEYWORD.findKey(k);
			if (v == FCS_KEYWORD.PnN) {
				paramNames.put(k, raw.get(k));
			}
			// System.out.println(k + "\t" + (v != null ? v.keyword : "null") + "\t" + raw.get(k));
		}

		setParameterNames(new String[paramNames.size()]);
		parameterNamesLong = new String[paramNames.size()];
		for (int i = 1; i < getParameterNames().length + 1; i++) {
			getParameterNames()[i - 1] = paramNames.get(FCS_KEYWORD.PnN.getKey(i));
			parameterNamesLong[i
												 - 1] = hasKeyword(FCS_KEYWORD.PnS.getKey(i)) ? getKeyword(FCS_KEYWORD.PnS.getKey(i))
																																			: getParameterNames()[i - 1];
		}

	}

	String[] getParameterNames() {
		return parameterNames;
	}

	void setParameterNames(String[] parameterNames) {
		this.parameterNames = parameterNames;
	}

}
