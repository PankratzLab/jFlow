package org.genvisis.one.ben.fcs.gating;

import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.genvisis.common.Logger;
import org.genvisis.common.Numbers;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class GateFileWriter {

	// built from: https://stackoverflow.com/questions/7373567/java-how-to-read-and-write-xml-files
	public static void writeGating(GatingStrategy gating, String outputFile, Logger log) {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder;
		try {
			builder = factory.newDocumentBuilder();
			Document doc = builder.newDocument();

			Element rootEle = doc.createElement("gating:Gating-ML");
			rootEle.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
			rootEle.setAttribute("xmlns:gating", "http://www.isac-net.org/std/Gating-ML/v2.0/gating");
			rootEle.setAttribute(	"xmlns:transforms",
														"http://www.isac-net.org/std/Gating-ML/v2.0/transformations");
			rootEle.setAttribute(	"xmlns:data-type",
														"http://www.isac-net.org/std/Gating-ML/v2.0/datatypes");
			rootEle.setAttribute(	"xsi:schemaLocation",
														"http://www.isac-net.org/std/Gating-ML/v2.0/gating http://www.isac-net.org/std/Gating-ML/v2.0/gating/Gating-ML.v2.0.xsd http://www.isac-net.org/std/Gating-ML/v2.0/transformations http://www.isac-net.org/std/Gating-ML/v2.0/gating/Transformations.v2.0.xsd http://www.isac-net.org/std/Gating-ML/v2.0/datatypes http://www.isac-net.org/std/Gating-ML/v2.0/gating/DataTypes.v2.0.xsd ");

			ArrayList<Gate> rootGates = gating.getRootGates();
			for (Gate g : rootGates) {
				addGateToDocument(g, doc, rootEle);
			}

			// Examples of opening tags for the two currently-supported types of gates
			// <gating:PolygonGate eventsInside="1" annoOffsetX="0" annoOffsetY="0" tint="#000000"
			// isTinted="0" lineWeight="Normal" userDefined="1" quadId="-1" gateResolution="256"
			// gating:id="ID793184926" >
			// <gating:RectangleGate eventsInside="1" annoOffsetX="0" annoOffsetY="0" tint="#000000"
			// isTinted="0" lineWeight="Normal" userDefined="1" percentX="0" percentY="0"
			// gating:id="ID1395151554" gating:parent_id="ID793184926" >

			doc.appendChild(rootEle);
			try {
				Transformer tr = TransformerFactory.newInstance().newTransformer();
				tr.setOutputProperty(OutputKeys.INDENT, "yes");
				tr.setOutputProperty(OutputKeys.METHOD, "xml");
				tr.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
				tr.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

				// send DOM to file
				tr.transform(new DOMSource(doc), new StreamResult(new FileOutputStream(outputFile)));

				log.report("Gating written to file: " + outputFile);
			} catch (TransformerException te) {
				log.reportException(te);
			} catch (IOException ioe) {
				log.reportException(ioe);
			}
		} catch (ParserConfigurationException e) {
			log.reportException(e);
		}


	}

	private static void addGateToDocument(Gate g, Document d, Element parentEle) {
		Element e;
		e = d.createElement(g.getXMLTag());

		if (g instanceof PolygonGate) {
			e.setAttribute("quadId", "-1"); // TODO "SpiderGates"
			e.setAttribute("gateResolution", "256"); // TODO custom info to note FlowJo inaccuracy?
		}

		e.setAttribute("gating:id", g.getID());
		if (g.getParentGate() != null) {
			e.setAttribute("gating:parent_id", g.getParentGate().getID());
		}

		for (GateDimension gd : g.getDimensions()) {
			Element dim1 = d.createElement("gating:dimension");

			if (gd instanceof RectangleGateDimension) {
				RectangleGateDimension rgd = (RectangleGateDimension) gd;
				if (Numbers.isFinite(rgd.getMax())) {
					dim1.setAttribute("gating:max", "" + Math.max(((RectangleGateDimension) gd).getMin(),
																												((RectangleGateDimension) gd).getMax()));
				}
				if (Numbers.isFinite(rgd.getMin())) {
					dim1.setAttribute("gating:min", "" + Math.min(((RectangleGateDimension) gd).getMin(),
																												((RectangleGateDimension) gd).getMax()));
				}
			}

			Element dim2 = d.createElement("data-type:fcs-dimension"); // TODO what constitutes a
																																	// "new-dimension"?
			dim2.setAttribute("data-type:name", gd.getParam());
			dim1.appendChild(dim2);
			e.appendChild(dim1);
		}

		if (g instanceof PolygonGate) {
			Path2D path = ((PolygonGate) g).getPath();
			PathIterator pi = path.getPathIterator(null);
			while (!pi.isDone()) {
				double[] coords = new double[6];
				int type = pi.currentSegment(coords);
				if (type == PathIterator.SEG_CLOSE) {
					pi.next();
					continue;
				}
				Element vert = d.createElement("gating:vertex");
				Element vert1 = d.createElement("gating:coordinate");
				Element vert2 = d.createElement("gating:coordinate");

				vert1.setAttribute("data-type:value", coords[0] + "");
				vert2.setAttribute("data-type:value", coords[1] + "");

				vert.appendChild(vert1);
				vert.appendChild(vert2);
				e.appendChild(vert);
				pi.next();
			}
		}

		parentEle.appendChild(e);

		for (Gate g2 : g.getChildGates()) {
			addGateToDocument(g2, d, parentEle);
		}

		// TODO these aren't necessary to reload in FlowJo, but perhaps should be included anyway?
		// eventsInside="1"
		// annoOffsetX="0"
		// annoOffsetY="0"
		// tint="#000000"
		// isTinted="0"
		// lineWeight="Normal"
		// userDefined="1"
		// percentX="0"
		// percentY="0"

	}

}
