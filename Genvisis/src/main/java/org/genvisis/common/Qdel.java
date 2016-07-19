package common;

public class Qdel {
	public static void main(String[] args) {
		int start, stop;
		
		if (args.length != 2) {
			System.err.println("Error - requires exactly 2 arguments (the first and last job to qdel)");
			return;
		}
		
		try {
			start = Integer.parseInt(args[0]);
        } catch (Exception e) {
	        System.err.println("Error - not a valid integer: "+args[0]);
	        return;
        }

        try {
			stop = Integer.parseInt(args[1]);
        } catch (Exception e) {
	        System.err.println("Error - not a valid integer: "+args[1]);
	        return;
        }
		
        for (int i = start; i<=stop; i++) {
            try {
    	        CmdLine.run("qdel "+i, "./");
            } catch (Exception e) {
    	        e.printStackTrace();
            }
        }

	}
}
