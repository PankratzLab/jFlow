package org.genvisis.widgets;

import java.util.Vector;

import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class ClassFromVariables {
  public static final boolean DEFAULT_NO_SET = true;

  public static void create() {
    String[] line, clip;
    String trans;
    Vector<String[]> v = new Vector<String[]>();
    String className = "Constructor";
    String[][] vars;
    boolean noset = DEFAULT_NO_SET;
    boolean serializable = false;

    clip = ext.getClipboard().split("\\n");
    trans = "";

    for (String element : clip) {
      line = ext.replaceAllWith(element, "{", "").trim().split("[\\s]+");
      if (line[0].equals("noset") || line[0].equals("no set")) {
        noset = true;
      } else if (line[0].equals("set")) {
        noset = false;
      } else if (element.contains("serialVersionUID")) {
        serializable = true;
      } else {
        trans += element + "\n";
        if (line.length == 1
            && (line[0].equals("") || line[0].equals("{") || line[0].equals("}"))) {

        } else if (line.length >= 3 && line[1].equals("class")) {
          className = line[2];
          if (line.length >= 5 && line[4].equals("Serializable")) {
            serializable = true;
            trans += "	public static final long serialVersionUID = 1L;\n\n";
          }
        } else if (line.length == 3) {
          if (!line[2].endsWith(";")) {
            System.err.println("Error - don't know how to parse: '" + element + "'");
          } else {
            line[2] = line[2].substring(0, line[2].length() - 1);
            v.add(line);
          }
        } else {
          System.err.println("Error - don't know how to parse: '" + element + "'");
        }
      }
    }

    vars = Matrix.toStringArrays(v);

    trans += "\n";
    trans += "\tpublic " + className + "(";
    for (int i = 0; i < vars.length; i++) {
      trans += (i == 0 ? "" : ", ") + vars[i][1] + " " + vars[i][2];
    }
    trans += ") {\n";
    for (String[] var : vars) {
      trans += "\t\tthis." + var[2] + " = " + var[2] + ";\n";
    }
    trans += "\t}\n";
    for (String[] var : vars) {
      if (!noset) {
        trans += "\n";
        trans += "\tpublic void set" + ext.capitalizeFirst(var[2]) + "(" + var[1] + " " + var[2]
                 + ") {\n";
        trans += "\t\tthis." + var[2] + " = " + var[2] + ";\n";
        trans += "\t}\n";
      }
      trans += "\n";
      trans += "\tpublic " + var[1] + " get" + ext.capitalizeFirst(var[2]) + "() {\n";
      trans += "\t\treturn " + var[2] + ";\n";
      trans += "\t}\n";
    }
    if (serializable) {
      trans += "\n";
      trans += "\tpublic void serialize(String filename) {\n";
      trans += "\t\tFiles.writeSerial(this, filename);\n";
      trans += "\t}\n";
      trans += "\n";
      trans += "\tpublic static " + className + " load(String filename, boolean jar) {\n";
      trans += "\t\treturn (" + className + ")Files.readSerial(filename, jar);\n";
      trans += "\t}";
    }

    ext.setClipboard(trans);
  }

  public static void main(String[] argv) {
    create();
  }

}
