package gwas;

public class Permute {


	public static void permute(String phenoFile) {
		
	}	
	
	
	
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String pheno = "pheno.dat";

		String usage = "\n" + "gwas.Permute requires 0-1 arguments\n"
				+ "   (1) filename (i.e. pheno=" + pheno + " (default))\n"
				+ "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				pheno = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			permute(pheno);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
