package widgets;

import common.*;

public class ClipSwap {
	public static final String CONTRACT_EXPAND_DELIM = "^";

	public static void contract() {
		String str = ext.getClipboard();
		
		if (CONTRACT_EXPAND_DELIM.length() == 1) {
			str = str.replace('\t', CONTRACT_EXPAND_DELIM.charAt(0));
		} else {
			str = ext.replaceAllWith(str, "\t", CONTRACT_EXPAND_DELIM);
		}
		
		ext.setClipboard(str);
	}

	public static void expand() {
		String str = ext.getClipboard();
		
		if (CONTRACT_EXPAND_DELIM.length() == 1) {
			str = str.replace(CONTRACT_EXPAND_DELIM.charAt(0), '\t');
		} else {
			str = ext.replaceAllWith(str, CONTRACT_EXPAND_DELIM, "\t");
		}
		
		ext.setClipboard(str);
	}
	
	public static void slash() {
		String str = ext.getClipboard();
		str = ext.replaceAllWith(str, "\\\\", "\\");
		str = ext.replaceAllWith(str, "\\", "qwerty");
//		str = ext.replaceAllWith(str, "qwerty", "\\\\");
//		ext.setClipboard(str+"\\\\");
		str = ext.replaceAllWith(str, "qwerty", "/");
		ext.setClipboard(str+"/");
	}

	public static void unique() {
		ext.setClipboard(Unique.proc(ext.getClipboard().trim().split("\\n")));
	}
	
	public static void removeFormatting() {
		ext.setClipboard(ext.getClipboard());
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    boolean slash = false;
	    boolean unique = false;
	    boolean contract = false;
	    boolean expand = false;
	    boolean removeFormatting = false;

	    String usage = "\n"+
	    "widgets.ClipSwap requires 0-1 arguments\n"+
	    "   (1) Fix slashes (i.e. -slash (not the default))\n"+
	    "   (2) Contracts contents of clipboard (i.e. -contract (not the default))\n"+
	    "   (3) Expands contents of clipboard (i.e. -expand (not the default))\n"+
	    "   (4) Find unique set and count counts (i.e. -unique (not the default))\n"+
	    "   (4) Remove formatting, leaving only plain text (i.e. -removeFormatting (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("-slash")) {
			    slash = true;
			    numArgs--;
		    } else if (args[i].startsWith("-contract")) {
			    contract = true;
			    numArgs--;
		    } else if (args[i].startsWith("-expand")) {
			    expand = true;
			    numArgs--;
		    } else if (args[i].startsWith("-unique")) {
			    unique = true;
			    numArgs--;
		    } else if (args[i].startsWith("-removeFormatting")) {
			    removeFormatting = true;
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
	    	if (slash) {
	    		slash();
	    	}
	    	if (contract) {
	    		contract();
	    	}
	    	if (expand) {
	    		expand();
	    	}
	    	if (unique) {
	    		unique();
	    	}
	    	if (removeFormatting) {
	    		removeFormatting();
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
