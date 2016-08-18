// calculates information content similar to mapmaker
package org.genvisis.link.init;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class info {

	public info(int start, int stop) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		PrintWriter markerWriter = null;
		PrintWriter masterWriter = new PrintWriter(new FileWriter("master info.out"));
		PrintWriter masterMarkerWriter = new PrintWriter(new FileWriter("master marker info.out"));
		StringTokenizer st;
		String chrome, trav, prev;
		Vector<String> positron = new Vector<String>();
		Vector<String> victor = new Vector<String>();
		Vector<String> zeroGoodness = new Vector<String>();
		Vector<String> markerNames = new Vector<String>();
		double pos, inf;

		for (int chromosome = start; chromosome<=stop; chromosome++) {
			if (chromosome<10) {
				chrome = "0"+chromosome;
			} else {
				chrome = ""+chromosome;
			}
			positron.removeAllElements();
			victor.removeAllElements();
			zeroGoodness.removeAllElements();
			markerNames.removeAllElements();
			try {
				try {
					reader = new BufferedReader(new FileReader("chromf"+chrome+".lin.out"));
				} catch (Exception e) {
					// try {
					// reader = new BufferedReader(new
					// FileReader("chromf"+chrome+".D.out"));
					// } catch (Exception e2) {
					// try {
					// reader = new BufferedReader(new
					// FileReader("chromf"+chrome+".R.out"));
					// } catch (Exception e3) {
					// System.err.println("Error: Not chromf"+chrome+".lin.out,
					// chromf"+chrome+".D.out or chromf"+chrome+".R.out could be
					// found for questioning.");
					// System.exit(0);
					// }
					// }
					System.err.println("Error: chromf"+chrome+".lin.out could be found for questioning.");
					System.exit(0);
				}

				reader.readLine();
				writer = new PrintWriter(new FileWriter("info"+chrome+".out"));
				markerWriter = new PrintWriter(new FileWriter("infoJustMarkers"+chrome+".out"));

				boolean first = true;
				int numFams = 1;
				int count = 0;
				prev = "";
				while (reader.ready()) {
					st = new StringTokenizer(reader.readLine());
					trav = st.nextToken();
					if (victor.isEmpty()) {
						prev = trav;
					}
					if (trav.equals(prev)) {
						count++;
						pos = Double.valueOf(st.nextToken()).doubleValue();
						st.nextToken();
						st.nextToken();
						st.nextToken();
						inf = Double.valueOf(st.nextToken()).doubleValue();
						if (first) {
							victor.add(inf+"");
							positron.add(pos+"");
							zeroGoodness.add("0");
							markerNames.add(st.nextToken());
						} else {
							victor.setElementAt((Double.valueOf(victor.elementAt(count)).doubleValue()+inf)+"", count);
						}
						// if (inf < 0.001) {
						// zeroGoodness.insertElementAt((Integer.valueOf(zeroGoodness.elementAt(count)).intValue()+1)+"",
						// count);
						// }
					} else {
						count = 0;
						numFams++;
						st.nextToken();
						st.nextToken();
						st.nextToken();
						st.nextToken();
						inf = Double.valueOf(st.nextToken()).doubleValue();
						// if (inf < 0.001) {
						// zeroGoodness.insertElementAt((Integer.valueOf(zeroGoodness.elementAt(count)).intValue()+1)+"",
						// count);
						// }
						victor.setElementAt((Double.valueOf(victor.elementAt(count)).doubleValue()+inf)+"", 0);
						first = false;
						prev = trav;
					}
				}
				reader.close();

				for (int i = 0; i<victor.size(); i++) {
					// writer.println(positron.elementAt(i)+"\t"+ext.formDeci(Double.valueOf(victor.elementAt(i)).doubleValue()/(numFams-Integer.valueOf(zeroGoodness.elementAt(count)).intValue()),
					// 5));
					writer.println(positron.elementAt(i)+"\t"+ext.formDeci(Double.valueOf(victor.elementAt(i)).doubleValue()/numFams, 5));
					masterWriter.println(positron.elementAt(i)+"\t"+ext.formDeci(Double.valueOf(victor.elementAt(i)).doubleValue()/numFams, 5));
					if (!markerNames.elementAt(i).equals("-")) {
						markerWriter.println(markerNames.elementAt(i)+"\t"+ext.formDeci(Double.valueOf(victor.elementAt(i)).doubleValue()/numFams, 5));
						masterMarkerWriter.println(markerNames.elementAt(i)+"\t"+ext.formDeci(Double.valueOf(victor.elementAt(i)).doubleValue()/numFams, 5));
					}
				}

			} catch (Exception e) {
				System.out.println("**crapped out on chromosome "+chromosome);
				// e.printStackTrace();
			}
			writer.close();
			markerWriter.close();
			System.out.println("Completed chromosome "+chromosome);
		}
		masterWriter.close();
		masterMarkerWriter.close();
	}

	public static void main(String[] args) throws IOException {
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-?")||args[i].equals("-help")||args[i].equals("-h")) {
				System.out.println("Expecting 0-2 arguments (chromosome number to start [and stop], default all).");
				System.exit(1);
			}
		}

		try {
			if (args.length==0) {
				new info(1, 23);
			}
			if (args.length==1) {
				new info(Integer.valueOf(args[0]).intValue(), Integer.valueOf(args[0]).intValue());
			}
			if (args.length==2) {
				new info(Integer.valueOf(args[0]).intValue(), Integer.valueOf(args[1]).intValue());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
