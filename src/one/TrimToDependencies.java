package one;

import java.io.*;
import java.util.*;
import common.*;

public class TrimToDependencies implements Serializable {
	private static final long serialVersionUID = 1L;

	public static final String DEPENDENCY_FINDER_BINARIES = "D:/home/npankrat/jProjects/DependencyFinder/bin/";
	public static final String[] BASE_CLASSES = {"Math"};

	private String[] filenames;
	private FileInfo[] fileBits;
	private Hashtable<String,Class> allClasses;
	private Hashtable<String,Method> allMethods;
	private Hashtable<String,Variable> allVariables;
	private Logger log;
	
	class FileInfo implements Serializable {
		private static final long serialVersionUID = 1L;
		private String name;
		private String packageName;
		private Vector<String> imports;
		private Vector<Class> classes;
		private boolean required;
		
		public FileInfo(String root, String filename) {
			BufferedReader reader;
			String[] line;
			String temp, trav;
			Class clas;
			boolean withinComment;
			String currentComment;
			boolean done;
			
			packageName = null;
			name = ext.replaceAllWith(filename.substring(0, filename.indexOf(".java")), new String[][] {{"/", "."}, {"\\", "."}});
			required = false;
			
			done = false;
			currentComment = null;
			withinComment = false;
			try {
				reader = new BufferedReader(new FileReader(root+filename));
				imports = new Vector<String>();
				classes = new Vector<TrimToDependencies.Class>();
				while (reader.ready()) {
					reader.mark(10000);
					temp = removeAnyTrailingComment(reader.readLine());
					if (temp.trim().startsWith("package ")) {
						line = temp.trim().split("[\\s]+");
						if (!line[1].substring(line[1].length()-1).equals(";")) {
							log.reportError("Error with package declaration");
						}
						packageName = line[1].substring(0, line[1].length()-1);
					} else if (temp.contains("import ")) {
						imports.add(temp);
					} else if (temp.contains("class ") || temp.contains("interface ")) {
						line = temp.trim().split("[\\s]+");
						reader.reset();
						trav = (packageName==null?"":packageName+".")+line[Math.max(ext.indexOfStr("class", line), ext.indexOfStr("interface", line)) +1];
						if (allClasses.containsKey(trav)) {
							clas = allClasses.get(trav);
						} else {
							allClasses.put(trav, clas = new Class(trav));
						}
						clas.populate(reader);
						clas.setParentFile(this);
						if (currentComment != null) {
							clas.setComment(currentComment);
							currentComment = null;
						}
						classes.add(clas);
					} else if (temp.trim().startsWith("//")){
					} else if (temp.contains("/*")){
						currentComment = temp+"\r\n";
						withinComment = true;
					} else if (withinComment) {
						currentComment += temp+"\r\n";
						if (temp.contains("*/")) {
							withinComment = false;
						}
					} else if (temp.contains("}")) {
						done = true;
					} else if (temp.trim().equals("")){
					} else if (done) {
						log.reportError("Error - line came after final bracket: "+temp);
					} else if (!reader.ready() && !done) {
						log.reportError("Error - never got final bracket in "+name);
					} else {
						log.reportError("Error - what to do with the following line from "+filename+":");
						log.reportError("\""+temp+"\"");
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + filename + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + filename + "\"");
				System.exit(2);
			}
		}
		
		public String getName() {
			return name;
		}

		public String getPackageName() {
			return packageName;
		}

		public Vector<String> getImports() {
			return imports;
		}

		public Vector<Class> getClasses() {
			return classes;
		}

		public void setRequired(boolean status) {
			required = status;
		}

		public boolean getRequired() {
			return required;
		}
	}
	
	class Variable implements Serializable {
		private static final long serialVersionUID = 1L;
		private String type;
		private String name;
		private boolean constant;
		private boolean required;
		
		public Variable(String newName) {
			name = newName;
			type = null;
			constant = false;
			required = false;
		}
		
		public void setAsConstant() {
			constant = true;
		}
		
		public void setType(String newType) {
			type = newType;
		}
		
		public String getType() {
			return type;
		}

		public String getName() {
			return name;
		}

		public boolean getConstant() {
			return constant;
		}

		public boolean getRequired() {
			return required;
		}
	}
	
	class Class implements Serializable {
		private static final long serialVersionUID = 1L;
		private String name;
		private String comment;
		private FileInfo parentFile;
		private Vector<Class> innerClasses;
		private Vector<Variable> variables;
		private Vector<Method> methods;
		private boolean required;
		
		public Class(String newName) {
			name = newName;
			parentFile = null;
			innerClasses = new Vector<Class>();
			methods = new Vector<Method>();
			variables = new Vector<Variable>();
			required = false;
		}
		
		public void populate(BufferedReader reader) {
			String[] line;
			String temp, trav;
			int numBrackets;
			Class clas;
			Method method;
			Variable variable;
			boolean withinComment;
			String currentComment;
			
//			, FileInfo fileInfo
//			if (parentFile == null) {
//				parentFile = fileInfo;
//			} else if (parentFile != fileInfo) {
//				System.err.println("Error - two different files are declaring a class called '"+name+"'");
//			}
			
			currentComment = null;
			withinComment = false;
			numBrackets = -9;
			try {
				log.report("Start of class '"+name+"'");
				temp = removeAnyTrailingComment(reader.readLine());
				if (temp.contains("extends")) {
					// oh well
				}
				numBrackets = updateBracketCount(numBrackets, temp);
				log.report(" C"+numBrackets+"\t"+temp);
				
				while (numBrackets == -9 || numBrackets > 0) {
					reader.mark(10000);
					temp = removeAnyTrailingComment(getRidOfAnythingInQuotes(reader.readLine()));
					line = temp.trim().split("[\\s]+");
					if (temp.contains("class ") || temp.contains("interface ")) {
						line = temp.trim().split("[\\s]+");
						reader.reset();
						trav = name+"."+line[Math.max(ext.indexOfStr("class", line), ext.indexOfStr("interface", line))+1];
						if (allClasses.containsKey(trav)) {
							clas = allClasses.get(trav);
						} else {
							allClasses.put(trav, clas = new Class(trav));
						}
						clas.populate(reader);
						clas.setParentFile(getParentFile());
						if (currentComment != null) {
							clas.setComment(currentComment);
							currentComment = null;
						}
						innerClasses.add(clas);
					} else if ((line[0].equals("public") || line[0].equals("private") || line[0].equals("protected")) && getRidOfAnythingInQuotes(temp).contains("(") && !getRidOfAnythingInQuotes(temp).contains("=")) {
						reader.reset();
						line = flipGenPars(ext.replaceAllWith(temp.substring(0, temp.indexOf("(")), new String[][] {{" static ", ""}, {"abstract ", ""}, {"final ", ""}, {"public ", ""}, {"private ", ""}, {"protected ", ""}})).trim().split("[\\s]+");
						if (line.length > 2) {
							log.reportError("Error - thought I was parsing a method here:");
							log.reportError(temp);
							System.exit(1);
						}
						trav = name+"."+line[line.length-1];
						if (allMethods.containsKey(trav)) {
							method = allMethods.get(trav);
						} else {
							allMethods.put(trav, method = new Method(trav));
						}
						if (line.length == 2) {
							method.setReturnType(line[0]);
						}
						method.populate(reader);
						if (currentComment != null) {
							method.setComment(currentComment);
							currentComment = null;
						}
						method.setParentClass(this);
						methods.add(method);
					} else if (temp.trim().startsWith("//")){
					} else if (temp.contains("/*")){
						currentComment = temp+"\r\n";
						withinComment = true;
						if (temp.contains("*/")) {
							withinComment = false;
						}
					} else if (withinComment) {
						currentComment += temp+"\r\n";
						if (temp.contains("*/")) {
							withinComment = false;
						}
					} else if (temp.contains("}") && updateBracketCount(0, temp) < 0) {
						numBrackets = updateBracketCount(numBrackets, temp);
					} else if (temp.trim().equals("")){
					} else if (temp.trim().equals("@Override") || temp.trim().startsWith("@SuppressWarnings")){
						currentComment += temp+"\r\n";
					} else {
						reader.reset();
						temp = getNextCompleteLine(reader, log);
						line = temp.trim().split("[\\s]+");
						if (!line[0].equals("public") && !line[0].equals("private") && !line[0].equals("protected")) {
							log.reportError("Error - assuming the following is an unsafe variable declaration; shame on you!");
							log.reportError(temp);
							System.exit(1);
						}
						line = ext.replaceAllWith(temp, new String[][] {{" static ", " "}, {"final ", ""}, {"public ", ""}, {"private ", ""}, {"protected ", ""}}).trim().split("[\\s]+");
						if (line.length < 2) {
							log.reportError("Error - assumed the following was a variable, but I guess not:");
							log.reportError(temp);
						}
						trav = name+"."+line[1];
						if (allVariables.containsKey(trav)) {
							variable = allVariables.get(trav);
						} else {
							allVariables.put(trav, variable = new Variable(trav));
						}
						variable.setType(line[0]);
						if (temp.contains("final ") || temp.contains("static ")) {
							variable.setAsConstant();
						}
						variables.add(variable);
						numBrackets = updateBracketCount(numBrackets, temp);
						log.report(" C"+numBrackets+"\t"+temp);
					}					
					
				}
				log.report("End of class '"+name+"'");
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + name + "\"");
				log.reportException(ioe);
				System.exit(2);
			}
			if (numBrackets != 0) {
				log.reportError("Error - must be some double brackets totalling "+numBrackets+" at the end of class '"+name+"'");
			}
		}
		
		public void determineDependents() {
			Method method;
//			Variable variable;
			
			for (int i = 0; i < methods.size(); i++) {
				method = methods.elementAt(i);
//				System.out.println(method.name);
				method.determineDeps();
			}

		}
		
		public String getName() {
			return name;
		}

		public void setComment(String currentComment) {
			comment = currentComment;
		}
		
		public String getComment() {
			return comment;
		}

		public void setParentFile(FileInfo fileInfo) {
			parentFile = fileInfo;
		}
		
		public FileInfo getParentFile() {
			return parentFile;
		}

		public void setRequired(boolean status) {
			required = status;
		}
		
		public boolean getRequired() {
			return required;
		}
		
	}

	class Method implements Serializable {
		private static final long serialVersionUID = 1L;
		private String name;
		private Class parentClass;
		private Hashtable<String,Variable> variables;
		private String returnType;
		private Vector<Method> methodsUsed;
		private String comment;
		private boolean required;
		private String[] allLines;
		
		public Method(String newName) {
			name = newName;
			variables = new Hashtable<String, Variable>();
			returnType = null;
			methodsUsed = new Vector<Method>();
			required = false;
			allLines = null;
		}
		
		public void populate(BufferedReader reader) {
			String[] line, sub;
			String temp;
			int numBrackets;
//			Method method;
			String varName;
			Variable variable;
			Vector<String> vLines;
			
			vLines = new Vector<String>();
			numBrackets = -9;
			try {
				log.report("Start of method '"+name+"'");
				temp = removeAnyTrailingComment(reader.readLine());
				numBrackets = updateBracketCount(numBrackets, temp);
				log.report("  M"+numBrackets+"\t"+temp);
				temp = temp.substring(temp.indexOf("(")+1, temp.indexOf(")"));
				if (temp.length() > 0) {
					line = flipGenPars(temp.trim()).split(",");
					for (int i = 0; i < line.length; i++) {
						sub = ext.replaceAllWith(line[i], new String[][] {{" static ", " "}, {"final ", ""}}).trim().split("[\\s]+");
						if (sub.length != 2) {
							log.reportError("Error - in argument: "+line[i]);
						}
						varName = sub[1];
						variables.put(varName, variable = new Variable(varName));
						variable.setType(sub[0]);
					}
				}

				while (numBrackets == -9 || numBrackets > 0) {
//					temp = getNextCompleteLine(reader);
					temp = removeAnyTrailingComment(getRidOfAnythingInQuotes(reader.readLine()));
//					line = temp.trim().split("[\\s]+");
//					if ((line[0].equals("public") || line[0].equals("private") || line[0].equals("protected")) && temp.contains("(")) {
//						reader.reset();
//						line = ext.replaceAllWith(temp, new String[][] {{" static ", ""}, {"(", " ("}}).trim().split("[\\s]+");
//						trav = line[1];
//						if (allMethods.containsKey(trav)) {
//							method = allMethods.get(trav);
//						} else {
//							allMethods.put(trav, method = new Method(trav));
//						}
//						method.populate(reader);
//						methods.add(method);
//					} else if (!temp.trim().equals("")){
//						System.err.println("Error - what to do with the following line from "+filename+":");
//						System.err.println("\""+temp+"\"");
//					}
//					
					vLines.add(temp);
					numBrackets = updateBracketCount(numBrackets, temp);
					log.report("  M"+numBrackets+"\t"+temp);
				}
				log.report("End of method '"+name+"'");
				allLines = Array.toStringArray(vLines);
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + name + "\"");
				log.reportException(ioe);
				System.exit(2);
			}
			if (numBrackets != 0) {
				log.reportError("Error - must be some double brackets totalling "+numBrackets+" at the end of class '"+name+"'");
			}
		}
		
		public void determineDeps() {
			String[] line;
			
			for (int i = 0; i < allLines.length; i++) {
				line = allLines[i].trim().split("[\\s]+");
				for (int j = 0; j < line.length; j++) {
					if (line[j].contains("(")) {
						if (j>0 && line[j-1].equals("new")) {
							System.out.println("\tConstructor: "+line[j]);
						} else {
							System.out.println("\t"+line[j]);
						}
					}
				}
			}
		}
		
		public String getName() {
			return name;
		}

		public void setComment(String currentComment) {
			comment = currentComment;
		}
		
		public String getComment() {
			return comment;
		}

		public void setParentClass(Class parent) {
			parentClass = parent;
		}
		
		public Class getParentClass() {
			return parentClass;
		}

		public void setReturnType(String type) {
			returnType = type;
		}

		public String getReturnType() {
			return returnType;
		}

		public Vector<Method> getMethodsUsed() {
			return methodsUsed;
		}

		public void setRequired(boolean status) {
			required = status;
		}
		
		public boolean getRequired() {
			return required;
		}
	}
	
	public static String removeAnyTrailingComment(String str) {
		if (str.contains("//")) {
			return str.substring(0, str.indexOf("//"));
		} else {
			return str;
		}
	}

	public static int updateBracketCount(int start, String str) {
		if (start == -9) {
			start = 0;
		}

		return start + ext.indicesWithinString("{", getRidOfAnythingInQuotes(str)).length - ext.indicesWithinString("}", getRidOfAnythingInQuotes(str)).length;
	}
	
	public static int netParenthesesCount(String str) {
		return ext.indicesWithinString("(", getRidOfAnythingInQuotes(str)).length - ext.indicesWithinString(")", getRidOfAnythingInQuotes(str)).length;
	}
	
	public static String flipGenPars(String str) {
		char[] chars;
		int depth;
		String newString;
		
		depth = 0;
		newString = "";
		chars = str.toCharArray();
		for (int i = 0; i < chars.length; i++) {
			if (depth > 0 && chars[i] == ' ') {
				
			} else {
				if (chars[i] == '<') {
					depth++;
				} else if (chars[i] == '>') {
					depth--;
				} else if (depth > 0 && chars[i] == ',') {
					chars[i] = '`';
				} else if (depth > 0 && chars[i] == '`') {
					chars[i] = ',';
				}
				newString += chars[i]+"";
			}
		}
		
		return newString;
	}
	
	public static String getRidOfAnythingInQuotes(String str) {
		char[] chars;
		boolean insideQuotes, insideQuote, fakeout;
		String newString;
		
		insideQuotes = false;
		insideQuote = false;
		fakeout = false;
		newString = "";
		chars = str.toCharArray();
		for (int i = 0; i < chars.length; i++) {
			if (fakeout) {
				fakeout = false;
			} else if (!insideQuote && chars[i] == '\"') {
				insideQuotes = !insideQuotes;
				newString += chars[i]+"";
			} else if (!insideQuotes && chars[i] == '\'') {
				insideQuote = !insideQuote;
				newString += chars[i]+"";
			} else if ((insideQuotes || insideQuote) && chars[i] == '\\') {
				fakeout = true;
			} else if (!insideQuotes && !insideQuote) {
				newString += chars[i]+"";
			}
		}
		
		return newString;
	}
	
	public static String getNextCompleteLine(BufferedReader reader, Logger log) {
		String str, temp;
		boolean acceptable;
		
		str = "";
		try {
			acceptable = false;
			do {
				temp = removeAnyTrailingComment(reader.readLine());
				str += temp+" ";
				acceptable = netParenthesesCount(str) == 0 && updateBracketCount(0, str) <= 0;
				if (netParenthesesCount(str) < 0) {
					log.reportError("Error - how can this be??");
					log.reportError(str);
					System.exit(1);
				}
			} while (!acceptable);
		} catch (IOException ioe) {
			log.reportError("Error reading file");
			System.exit(2);
		}
		
		return ext.replaceAllWith(str, new String[][] {{"\t\t", "\t"}, {"\t", " "}});
	}

	public TrimToDependencies(String source_directory, String target_directory) {
		allClasses = new Hashtable<String,Class>();
		allMethods = new Hashtable<String,Method>();
		allVariables = new Hashtable<String,Variable>();

		log = new Logger();
//		log = new Logger("trimming.log", false);

		filenames = Files.listAllFilesInTree(source_directory, false);
		fileBits = new FileInfo[filenames.length];
		for (int i = 0; i < filenames.length; i++) {
			if (filenames[i].endsWith(".java")) {
				fileBits[i] = new FileInfo(source_directory, filenames[i]);
			}
		}
	}
	
	public void linkTo(String[] coreClasses, String freshJarFile) {
		Class clas;
		
		if (new File(freshJarFile+".dependencies").exists()) {
			
			
			
			
			for (int i = 0; i < coreClasses.length; i++) {
				if (!allClasses.containsKey(coreClasses[i])) {
					System.err.println("Error - '"+coreClasses[i]+"' is not present in this massive archive");
					System.exit(1);
				}
				clas = allClasses.get(coreClasses[i]);
				System.out.println(clas.getName());
				clas.determineDependents();
			}
			
			
			
			
		} else {
			if (new File(freshJarFile).exists()) {
				System.err.println("The file '"+freshJarFile+".dependencies' has not yet been created; creating a batch file to process '"+freshJarFile+"' called "+freshJarFile+".genDeps.bat in the directory: "+ext.parseDirectoryOfFile(freshJarFile));
				if (ext.getTimeSince(new File(freshJarFile).lastModified(), 'H') < 3) {
					System.err.println("Warning - file '"+freshJarFile+"' is more than 3 hours old. Make sure it matches the source code exactly!");
				}
				Files.write("DependencyExtractor.bat "+freshJarFile+" -xml -maximize 1> "+freshJarFile+".dependencies", freshJarFile+".genDeps.bat");
			} else {
				System.err.println("The file '"+freshJarFile+".dependencies' has not yet been created, and the file '"+freshJarFile+"' could not be found in order to create it.");
			}
			
		}
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}
	
	public static TrimToDependencies load(String filename) {
		return (TrimToDependencies)Files.readSerial(filename, false, new Logger(), false);
	}

	public static void main(String[] args) {
		String source_dir = "D:/home/npankrat/jProjects/master/src/";
		String target_dir = "D:/home/npankrat/jProjects/target/src/";
		String[] coreClasses = new String[] {"filesys.DosageData"};
		String freshJarFile = "D:/home/npankrat/park.jar";
		TrimToDependencies trimmer;		

		try {
			if (new File("master.deps.ser").exists() && ext.getTimeSince(new File("master.deps.ser").lastModified(), 'D') < 7) {
				System.err.println("Warning - loading source code from memory. Delete/rename file 'master.deps.ser' if you want to recreate from scratch.");
				trimmer = load("master.deps.ser");
			} else {
				trimmer = new TrimToDependencies(source_dir, target_dir);
				trimmer.serialize("master.deps.ser");
			}
			trimmer.linkTo(coreClasses, freshJarFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
