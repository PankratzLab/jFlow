package org.genvisis.common;

import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;

public class HtmlSelection implements Transferable {
	private static ArrayList<DataFlavor> htmlFlavors;
	private final String html;

	public HtmlSelection(String html) {
		htmlFlavors = new ArrayList<DataFlavor>();
		try {
			htmlFlavors.add(new DataFlavor("text/html;class=java.lang.String"));
			htmlFlavors.add(new DataFlavor("text/html;class=java.io.Reader"));
			htmlFlavors.add(new DataFlavor("text/html;charset=unicode;class=java.io.InputStream"));
		} catch (ClassNotFoundException ex) {
			ex.printStackTrace();
		}

		this.html = html;
	}

	@Override
	public DataFlavor[] getTransferDataFlavors() {
		return htmlFlavors.toArray(new DataFlavor[htmlFlavors.size()]);
	}

	@Override
	public boolean isDataFlavorSupported(DataFlavor flavor) {
		return htmlFlavors.contains(flavor);
	}

	@Override
	public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException {
		if (String.class.equals(flavor.getRepresentationClass())) {
			return html;
		} else if (Reader.class.equals(flavor.getRepresentationClass())) {
			return new StringReader(html);
		}
		throw new UnsupportedFlavorException(flavor);
	}
}
