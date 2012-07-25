package cnv.plots;

public class deNovoCalling {
	byte[] motherCnvStatus, fatherCnvStatus, childCnvStatus;
	
	public deNovoCalling() {
		byte[] result;

		loadData();
		
		result = new byte[childCnvStatus.length];
		for (int i=0; i<childCnvStatus.length; i++) {
			if (motherCnvStatus[i]>0 && fatherCnvStatus[i]==0 && motherCnvStatus[i]==0) {
				result[i]=1;
			}
		}
		
	}
	
	public void loadData() {
		//TODO read from file
		motherCnvStatus=null;
		fatherCnvStatus=null;
		childCnvStatus=null;
		
		//TODO link the data by family ID
	}

}
